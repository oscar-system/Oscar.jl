using JSON
using UUIDs

include("serializers.jl")

const backref_sym = Symbol("#backref")
################################################################################
# Meta Data

@Base.kwdef struct MetaData
  author_orcid::Union{String, Nothing} = nothing
  name::Union{String, Nothing} = nothing
end

function metadata(;args...)
  return MetaData(;args...)
end

################################################################################
# Version info

function get_version_info()
  result = Dict{Symbol, Any}(
    :Oscar => ["https://github.com/oscar-system/Oscar.jl", VERSION_NUMBER]
  )
  return result
end
const oscarSerializationVersion = get_version_info()

################################################################################
# (De|En)coding types
#
# parameters of type should not matter here
const reverseTypeMap = Dict{String, Type}()

function registerSerializationType(@nospecialize(T::Type), str::String)
  if haskey(reverseTypeMap, str) && reverseTypeMap[str] != T
    error("encoded type $str already registered for a different type: $T versus $(reverseTypeMap[str])")
  end
  reverseTypeMap[str] = T
end

# registerSerializationType is a macro to ensure that the string we generate
# matches exactly the expression passed as first argument, and does not change
# in unexpected ways when import/export statements are adjusted.
# It also sets the value of serialize_with_id, which determines
# whether or not the type can be back referenced.
# If omitted, the default is that no back references are allowed.
function registerSerializationType(ex::Any,
                                   uses_id::Bool,
                                   str::Union{String,Nothing} = nothing)
  if str === nothing
    str = string(ex)
  end
  return esc(
    quote
      registerSerializationType($ex, $str)
      encode_type(::Type{<:$ex}) = $str
      # There exist types where equality cannot be discerned from the serialization
      # these types require an id so that equalities can be forced upon load.
      # The ids are only necessary for parent types, checking for element type equality
      # can be done once the parents are known to be equal.
      # For example two serializations of QQ[x] require ids to check for equality.
      # Although they're isomorphic rings, they may want to be treated as seperate
      # This is done since other software might not use symbols in their serialization of QQ[x].
      # Which will then still allow for the distinction between QQ[x] and QQ[y], i.e.
      # whenever there is a possibility (amongst any software system) that the objects
      # cannot be distinguish on a syntatic level we use ids.
      # Types like ZZ, QQ, and ZZ/nZZ do not require ids since there is no syntatic
      # ambiguities in their encodings.

      serialize_with_id(obj::T) where T <: $ex = $uses_id 
      serialize_with_id(T::Type{<:$ex}) = $uses_id 
    end)
end

macro registerSerializationType(ex::Any, str::Union{String,Nothing} = nothing)
  return registerSerializationType(ex, false, str)
end

macro registerSerializationType(ex::Any, uses_id::Bool, str::Union{String,Nothing} = nothing)
  return registerSerializationType(ex, uses_id, str)
end


function encode_type(::Type{T}) where T
  error("unsupported type '$T' for encoding")
end

function decode_type(input::String)
  get(reverseTypeMap, input) do
    error("unsupported type '$input' for decoding")
  end
end

function decode_type(input::Dict{Symbol, Any})
  return decode_type(input[:name])
end

################################################################################
# Encoding helper functions

@doc raw"""
    is_basic_serialization_type(::Type)

During the serialization of types of the form `Vector{T}`, entries
of type `T` will either be serialized as strings if `is_basic_serialization_type`
returns `true`, or serialized as a dict provided the serialization for such a `T`
exists. If `Vector{T}` is serialized with `is_basic_serialization_type(T) = true`
then the `entry_type` keyword is used to store the type `T` as a property of
the vector.

# Examples

```jldoctest
julia> is_basic_serialization_type(ZZRingElem)
true
```
"""
is_basic_serialization_type(::Type) = false
is_basic_serialization_type(::Type{ZZRingElem}) = true
is_basic_serialization_type(::Type{QQFieldElem}) = true
is_basic_serialization_type(::Type{Bool}) = true
is_basic_serialization_type(::Type{String}) = true
is_basic_serialization_type(::Type{Symbol}) = true
# this deals with int32, int64 etc.:precision
is_basic_serialization_type(::Type{T}) where T <: Number = isconcretetype(T)

# ATTENTION
# We need to distinguish between data with a globally defined normal form and data where such a normal form depends on some parameters.
# In particular, this does NOT ONLY depend on the type; see, e.g., FqField.

# Objects having a basic encoding are fundamental objects. They are "atomic"
# with regards to the available types in Oscar (maybe in general), they are base cases when
# a parent tree is constructed via recursion, and they can be decoded from a string or array of
# strings provided their parent is known. All types with basic serialization have a basic encoding
# has_elem_basic_encoding also deals with basic parametrized types.

# a basic serialization type is a type that only requires a single property for serialization.
# examples: QQ, Symbol, Int

# a basic encoded type is a type that only requires a single property or a list of ordered properties
# of the same type and a known parent for serialization.
# example: elements of ZZ/nZZ. Say the parent is ZZ/7ZZ then we can (de)serialize "6", in a sense once the parent is known
# the serialization relies on a "basic" serialization or list of ordered "basic" serializations. We use the word encoding
# here to mean serialization knowing how elements are represented in the parent.

# has_elem_basic_encoding is used for parent types (and only makes sense for parent types)

function has_elem_basic_encoding(obj::T) where T <: Ring
  return is_basic_serialization_type(elem_type(obj))
end

has_basic_encoding(obj::T) where T = is_basic_serialization_type(T)

function has_basic_encoding(obj::T) where T <: RingElem
  return has_elem_basic_encoding(parent(obj))
end

has_elem_basic_encoding(obj::T) where T = false

# used for types that require parents when serialized
type_needs_params(::Type) where Type = false

################################################################################
# High level

function load_ref(s::DeserializerState, id::String)
  if haskey(global_serializer_state.id_to_obj, UUID(id))
    loaded_ref = global_serializer_state.id_to_obj[UUID(id)]
  else
    ref_dict = s.refs[Symbol(id)]
    ref_dict[:id] = id
    loaded_ref = load_typed_object(s, ref_dict)
    global_serializer_state.id_to_obj[UUID(id)] = loaded_ref
  end
  return loaded_ref
end

function save_as_ref(s::SerializerState, obj::T) where T
  # find ref or create one
  ref = get(global_serializer_state.obj_to_id, obj, nothing)
  if ref !== nothing
    if !(ref in s.refs)
      push!(s.refs, ref)
    end
    return string(ref)
  end
  ref = global_serializer_state.obj_to_id[obj] = uuid4()
  global_serializer_state.id_to_obj[ref] = obj
  push!(s.refs, ref)
  return string(ref)
end

function save_object(s::SerializerState, x::Any, key::Symbol)
  s.key = key
  save_object(s, x)
end

function save_json(s::SerializerState, x::Any)
  data_json(s, x)
end

function save_json(s::SerializerState, x::Any, key::Symbol)
  s.key = key
  save_json(s, x)
end

function save_header(s::SerializerState, h::Dict{Symbol, Any}, key::Symbol)
  s.key = key
  data_dict(s) do
    for (k, v) in h
      save_object(s, v, k)
    end
  end
end

function save_typed_object(s::SerializerState, x::T) where T
  if type_needs_params(T)
    save_type_params(s, x, :type)
    save_object(s, x, :data)
  elseif Base.issingletontype(T)
    save_object(s, encode_type(T), :type)
  else
    save_object(s, encode_type(T), :type)
    save_object(s, x, :data)
  end
end

function save_typed_object(s::SerializerState, x::T, key::Symbol) where T
  s.key = key
  if serialize_with_id(x)
    # key should already be set before function call
    ref = save_as_ref(s, x)
    save_object(s, ref)
  else
    data_dict(s) do 
      save_typed_object(s, x)
    end
  end
end

function save_type_params(s::SerializerState, obj::Any, key::Symbol)
  s.key = key
  save_type_params(s, obj)
end

# The load mechanism first checks if the type needs to load necessary
# parameters before loading it's data, if so a type tree is traversed
function load_typed_object(s::DeserializerState, dict::Dict{Symbol, Any};
                           override_params::Any = nothing)
  T = decode_type(dict[:type])
  if Base.issingletontype(T) && return T()
  elseif type_needs_params(T)
    if !isnothing(override_params)
      params = load_type_params(s, T, override_params)
    else
      # depending on the type, :params is either an object to be loaded or a
      # dict with keys and object values to be loaded
      params = load_type_params(s, T, dict[:type][:params])
    end
    return load_object(s, T, dict[:data], params)
  else
    return load_object(s, T, dict[:data])
  end
end

function load_typed_object(s::DeserializerState, id::String)
  return load_ref(s, id)
end

################################################################################
# Default generic save_internal, load_internal
function save_object_generic(s::SerializerState, obj::T) where T
  s.key = :data
  data_dict(s) do
    for n in fieldnames(T)
      if n != :__attrs
        save_typed_object(s, getfield(obj, n), Symbol(n))
      end
    end
  end
end

function load_object_generic(s::DeserializerState, ::Type{T}, dict::Dict) where T
  fields = []
  for (n,t) in zip(fieldnames(T), fieldtypes(T))
    if n!= :__attrs
      println(t)
      push!(fields, load_object(s, t, dict[n]))
    end
  end
  return T(fields...)
end

################################################################################
# Utility functions for parent tree

# loads parent tree
function load_parents(s::DeserializerState, parent_ids::Vector)
  loaded_parents = []

  for id in parent_ids
    loaded_parent = load_ref(s, id)
    push!(loaded_parents, loaded_parent)
  end
  return loaded_parents
end

################################################################################
# Include serialization implementations for various types

include("basic_types.jl")
include("containers.jl")
include("PolyhedralGeometry.jl")
include("Combinatorics.jl")
include("Fields.jl")
include("ToricGeometry.jl")
include("Rings.jl")
include("Algebras.jl")
include("polymake.jl")
include("TropicalGeometry.jl")
include("QuadForm.jl")

################################################################################
# Include upgrade scripts

include("upgrades/main.jl")

function get_file_version(dict::Dict{Symbol, Any})
  ns = dict[:_ns]
  version_dict = ns[:Oscar][2]
  return get_version_number(version_dict)
end

function get_version_number(v_number::String)
  return VersionNumber(v_number)
end

################################################################################
# Interacting with IO streams and files

"""
    save(io::IO, obj::Any; metadata::MetaData=nothing)
    save(filename::String, obj::Any, metadata::MetaData=nothing)

Save an object `T` to the given io stream
respectively to the file `filename`.

See [`load`](@ref).

# Examples

```jldoctest
julia> meta = metadata(author_orcid="0000-0000-0000-0042", name="the meaning of life the universe and everything")
Oscar.MetaData("0000-0000-0000-0042", "the meaning of life the universe and everything")

julia> save("/tmp/fourtitwo.json", 42; metadata=meta);

julia> load("/tmp/fourtitwo.json")
42
```
"""
function save(io::IO, obj::T; metadata::Union{MetaData, Nothing}=nothing) where T
  state = serializer_open(io)
  data_dict(state) do
    save_typed_object(state, obj)

    if serialize_with_id(T)
      ref = get(global_serializer_state.obj_to_id, obj, nothing)
      if isnothing(ref)
        ref = global_serializer_state.obj_to_id[obj] = uuid4()
        global_serializer_state.id_to_obj[ref] = obj
      end
      save_object(state, string(ref), :id)
    end
    
    # this should be handled by serializers in a later commit / PR
    state.key = :refs
    !isempty(state.refs) && data_dict(state) do
      for id in state.refs
        ref_obj = global_serializer_state.id_to_obj[id]
        state.key = Symbol(id)
        data_dict(state) do
          save_typed_object(state, ref_obj)
        end
      end
    end
    save_header(state, oscarSerializationVersion, :_ns)
    
    if !isnothing(metadata)
      save_json(state, json(metadata), :meta)
    end
  end
  serializer_close(state)
  return nothing
end

function save(filename::String, obj::Any; metadata::Union{MetaData, Nothing}=nothing)
  temp_file = tempname(dirname(filename))
  open(temp_file, "w") do file
    save(file, obj; metadata=metadata)
  end
  Base.Filesystem.rename(temp_file, filename) # atomic "multi process safe"
  return nothing
end

"""
    load(io::IO; params::Any = nothing, type::Any = nothing)
    load(filename::String; params::Any = nothing, type::Any = nothing)

Load the object stored in the given io stream
respectively in the file `filename`.

If `params` is specified, then the root object of the loaded data
either will attempt a load using these parameters. In the case of Rings this
results in setting its parent, or in the case of a container of ring types such as
`Vector` or `Tuple`, then the parent of the entries will be set using their
 `params`.

If a type `T` is given then attempt to load the root object of the data
being loaded with this type; if this fails, an error is thrown.

See [`save`](@ref).

# Examples

```jldoctest
julia> save("/tmp/fourtitwo.json", 42);

julia> load("/tmp/fourtitwo.json")
42

julia> load("/tmp/fourtitwo.json"; type=Int64)
42

julia> R, x = QQ["x"]
(Univariate polynomial ring in x over QQ, x)

julia> p = x^2 - x + 1
x^2 - x + 1

julia> save("/tmp/p.json", p)

julia> p_loaded = load("/tmp/p.json", params=R)
x^2 - x + 1

julia> parent(p_loaded) === R
true

julia> save("/tmp/p_v.json", [p, p])

julia> loaded_p_v = load("/tmp/p_v.json", params=R)
2-element Vector{QQPolyRingElem}:
 x^2 - x + 1
 x^2 - x + 1

julia> parent(loaded_p_v[1]) === parent(loaded_p_v[2]) === R
true
```
"""
function load(io::IO; params::Any = nothing, type::Any = nothing)
  state = deserializer_open(io)

  # this should be moved to the serializer at some point
  jsondict = JSON.parse(io, dicttype=Dict{Symbol, Any})

  # handle different namespaces
  @req haskey(jsondict, :_ns) "Namespace is missing"
  _ns = jsondict[:_ns]
  if haskey(_ns, :polymake)
    # If this is a polymake file
    return load_from_polymake(jsondict)
  end
  @req haskey(_ns, :Oscar) "Not an Oscar object"

  # deal with upgrades
  file_version = get_file_version(jsondict)

  if file_version < VERSION_NUMBER
    jsondict = upgrade(jsondict, file_version)
  end

  # add refs to state for referencing during recursion
  if haskey(jsondict, :refs)
    merge!(state.refs, jsondict[:refs])
  end

  if type !== nothing
    # Decode the stored type, and compare it to the type `T` supplied by the caller.
    # If they are identical, just proceed. If not, then we assume that either
    # `T` is concrete, in which case `T <: U` should hold; or else `U` is
    # concrete, and `U <: T` should hold.
    #
    # This check should maybe change to a check on the whole type tree?
    U = decode_type(jsondict[:type])
    U <: type || U >: type || error("Type in file doesn't match target type: $(dict[:type]) not a subtype of $T")

    if type_needs_params(type)
      if isnothing(params) 
        params = load_type_params(state, type, jsondict[:type][:params])
      end
      return load_object(state, type, jsondict[:data], params)
    end
    Base.issingletontype(type) && return type()
    return load_object(state, type, jsondict[:data])
  end
  return load_typed_object(state, jsondict; override_params=params)
end

function load(filename::String; params::Any = nothing, type::Any = nothing)
  open(filename) do file
    return load(file; params=params, type=type)
  end
end
