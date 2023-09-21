################################################################################
# Description of the saving and loading mechanisms
#
# We require that any types serialized through OSCAR are registered using the
# @register_serialization_type macro. For more information, see its docstring.

# There are three pairs of saving and loading functions that are used
# during serialization:
# 1. save_typed_object, load_typed_object;
# 2. save_object, load_object.
# 3. save_type_params, load_type_params;
#
# In the following, we discuss each pair in turn.

################################################################################
# save_type_object / load_type_object

# For the most part these functions should not be touched, they are high level
# functions and are used to (de)serialize the object with its
# type information as well as its data. The data and type nodes are
# set in save_typed_object resulting in a "data branch" and "type branch".
# The usage of these functions can be used inside save_object / load_object
# and save_type_params / load_type_params. However using save_typed_object inside
# a save_object implementation will lead to a verbose format and should at some
# point be move to save_type_params.

################################################################################
# save_object / load_object

# These functions are at the core of the serialization and are the first functions
# that should be implemented when working on the serialization of a new type.
# Here is where one should use functions save_data_dict and save_data_array to structure
# the serialization. The examples show they can be used to save data using the structure
# of an array or dict. Each nested call to save_data_dict or save_data_array should be
# called with a key that can be passed as the second parameter.

#  Examples
#  function save_object(s::SerializerState, obj::NewType)
#    save_data_array(s) do
#      save_object(s, obj.1)
#      save_object(s, obj.2)
#
#      save_data_dict(s) do
#        save_object(s, obj.3, :key1)
#        save_object(s, obj.4, :key2)
#      end
#    end
#  end
#
#  This will result in a data format that looks like
#  [
#    obj.1,
#    obj.2,
#    {
#      "key1": obj.3,
#      "key2": obj.4
#    }
#  ]

#  function save_object(s::SerializerState, obj::NewType)
#    save_data_dict(s) do
#      save_object(s, obj.1, :key1)
#      save_data_array(s, :key2) do
#        save_object(s, obj.3)
#        save_typed_object(s, obj.4) # This is ok
#      end
#    end
#  end

#  This will result in a data format that looks like
#  {
#    "key1": obj.1,
#    "key2":[
#      obj.3,
#      {
#        "type": "Type of obj.4",
#        "data": obj.4
#      }
#    ]
#  }
#
# This is ok
# function save_object(s::SerializerState, obj:NewType)
#   save_object(s, obj.1)
# end

# This will throw an error 
# function save_object(s::SerializerState, obj:NewType)
#   save_object(s, obj.1, :key)
# end

# If you insist on having a key you should, first open a save_data_dict.

# function save_object(s::SerializerState, obj:NewType)
#   save_data_dict(s) do
#     save_object(s, obj.1, :key)
#   end
# end
#

# note for now save_typed_object must be wrapped in either a save_data_array or
# save_data_dict. Otherwise you will get a key override error.

# function save_object(s::SerializerState, obj:NewType)
#   save_data_dict(s) do
#     save_typed_object(s, obj.1, :key)
#   end
# end
#


################################################################################
# save_type_params / load_type_params

# The serialization mechanism stores data in the format of a tree, with the
# exception that some nodes may point to a shared reference. The "data branch"
# is anything that is a child node of a data node, whereas the "type branch" is
# any information that is stored in a node that is a child of a type node.
# Avoiding type information inside the data branch will lead to a more
# efficient serialization format. When the uses_params is set while calling
# @register_serialization_type (de)serialization will use 
# save_type_params / load_type_params to format the type information.
# In general we expect that implementing a save_type_params and load_type_params
# should not always be necessary. Many types will serialize their types
# in a similar fashion for example serialization of a FieldElem will
# use the save_type_params / load_type_params from RingElem since in both cases
# the only parameter needed for such types is their parent.

using JSON
using UUIDs

include("serializers.jl")

const type_key = :_type
const refs_key = :_refs
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
const oscar_serialization_version = get_version_info()

################################################################################
# (De|En)coding types

# parameters of type should not matter here
const reverse_type_map = Dict{String, Type}()

function register_serialization_type(@nospecialize(T::Type), str::String)
  if haskey(reverse_type_map, str) && reverse_type_map[str] != T
    error("encoded type $str already registered for a different type: $T versus $(reverse_type_map[str])")
  end
  reverse_type_map[str] = T
end

# @register_serialization_type NewType "String Representation of type" uses_id uses_params

# register_serialization_type is a macro to ensure that the string we generate
# matches exactly the expression passed as first argument, and does not change
# in unexpected ways when import/export statements are adjusted.
# The last three arguments are optional and can arise in any order. Passing a string
# argument will override how the type is stored as a string. The last two are boolean
# flags. When setting uses_id the object will be stored as a reference and will be
# referred to throughout the serialization using a UUID. This should typically only
# be used for types that do not have a fixed normal form for example PolyRing and MPolyRing.
# Using the uses_params flag will serialize the object with a more structured type
# description which will make the serialization more efficient see the discussion on
# save_type_params / load_type_params below.
function register_serialization_type(ex::Any, str::String, uses_id::Bool, uses_params::Bool)
  return esc(
    quote
      register_serialization_type($ex, $str)
      encode_type(::Type{<:$ex}) = $str
      # There exist types where equality cannot be discerned from the serialization
      # these types require an id so that equalities can be forced upon load.
      # The ids are only necessary for parent types, checking for element type equality
      # can be done once the parents are known to be equal.
      # For example two serializations of QQ[x] require ids to check for equality.
      # Although they're isomorphic rings, they may want to be treated as separate
      # This is done since other software might not use symbols in their serialization of QQ[x].
      # Which will then still allow for the distinction between QQ[x] and QQ[y], i.e.
      # whenever there is a possibility (amongst any software system) that the objects
      # cannot be distinguish on a syntactic level we use ids.
      # Types like ZZ, QQ, and ZZ/nZZ do not require ids since there is no syntactic
      # ambiguities in their encodings.

      serialize_with_id(obj::T) where T <: $ex = $uses_id 
      serialize_with_id(T::Type{<:$ex}) = $uses_id
      serialize_with_params(T::Type{<:$ex}) = $uses_params
    end)
end

macro register_serialization_type(ex::Any, args...)
  uses_id = false
  uses_params = false
  str = nothing
  for el in args
    if el isa String
      str = el
    elseif el == :uses_id
      uses_id = true
    elseif el == :uses_params
      uses_params = true
    end
  end
  if str === nothing
    str = string(ex)
  end

  return register_serialization_type(ex, str, uses_id, uses_params)
end

function encode_type(::Type{T}) where T
  error("unsupported type '$T' for encoding")
end

function decode_type(input::String)
  get(reverse_type_map, input) do
    error("unsupported type '$input' for decoding")
  end
end

function decode_type(input::Dict{Symbol, Any})
  return decode_type(input[:name])
end

# ATTENTION
# We need to distinguish between data with a globally defined normal form and data where such a normal form depends on some parameters.
# In particular, this does NOT ONLY depend on the type; see, e.g., FqField.

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
  set_key(s, key)
  save_object(s, x)
end

function save_json(s::SerializerState, x::Any)
  save_data_json(s, x)
end

function save_json(s::SerializerState, x::Any, key::Symbol)
  set_key(s, key)
  save_json(s, x)
end

function save_header(s::SerializerState, h::Dict{Symbol, Any}, key::Symbol)
  save_data_dict(s, key) do
    for (k, v) in h
      save_object(s, v, k)
    end
  end
end

function save_typed_object(s::SerializerState, x::T) where T
  if serialize_with_params(T)
    save_type_params(s, x, type_key)
    save_object(s, x, :data)
  elseif Base.issingletontype(T)
    save_object(s, encode_type(T), type_key)
  else
    save_object(s, encode_type(T), type_key)
    save_object(s, x, :data)
  end
end

function save_typed_object(s::SerializerState, x::T, key::Symbol) where T
  set_key(s, key)
  if serialize_with_id(x)
    # key should already be set before function call
    ref = save_as_ref(s, x)
    save_object(s, ref)
  else
    save_data_dict(s) do 
      save_typed_object(s, x)
    end
  end
end

function save_type_params(s::SerializerState, obj::Any, key::Symbol)
  set_key(s, key)
  save_type_params(s, obj)
end

# general loading of a reference
function load_type_params(s::DeserializerState, ::Type, ref::String)
  return load_ref(s, ref)
end

# The load mechanism first checks if the type needs to load necessary
# parameters before loading it's data, if so a type tree is traversed
function load_typed_object(s::DeserializerState, dict::Dict{Symbol, Any};
                           override_params::Any = nothing)
  T = decode_type(dict[type_key])
  if Base.issingletontype(T) && return T()
  elseif serialize_with_params(T)
    if !isnothing(override_params)
      params = override_params
    else
      # depending on the type, :params is either an object to be loaded or a
      # dict with keys and object values to be loaded
      params = load_type_params(s, T, dict[type_key][:params])
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
  save_data_dict(s, :data) do
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
  version_info = ns[:Oscar][2]
  return get_version_number(version_info)
end

function get_version_number(v_number::String)
  return VersionNumber(v_number)
end

# needed for older versions 
function get_version_number(dict::Dict)
  return VersionNumber(dict[:major], dict[:minor], dict[:patch])
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
  save_data_dict(state) do
    # write out the namespace first
    save_header(state, oscar_serialization_version, :_ns)

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
    !isempty(state.refs) && save_data_dict(state, refs_key) do
      for id in state.refs
        ref_obj = global_serializer_state.id_to_obj[id]
        state.key = Symbol(id)
        save_data_dict(state) do
          save_typed_object(state, ref_obj)
        end
      end
    end
    
    if !isnothing(metadata)
      save_json(state, json(metadata), :meta)
    end
  end
  serializer_close(state)
  return nothing
end

function save(filename::String, obj::Any; metadata::Union{MetaData, Nothing}=nothing)
  dir_name = dirname(filename)
  # julia dirname does not return "." for plain filenames without any slashes
  temp_file = tempname(isempty(dir_name) ? pwd() : dir_name)
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

  if haskey(jsondict, :id)
    id = jsondict[:id]
    if haskey(global_serializer_state.id_to_obj, UUID(id))
      return global_serializer_state.id_to_obj[UUID(id)]
    end
  end

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
  if haskey(jsondict, refs_key)
    merge!(state.refs, jsondict[refs_key])
  end

  if type !== nothing
    # Decode the stored type, and compare it to the type `T` supplied by the caller.
    # If they are identical, just proceed. If not, then we assume that either
    # `T` is concrete, in which case `T <: U` should hold; or else `U` is
    # concrete, and `U <: T` should hold.
    #
    # This check should maybe change to a check on the whole type tree?
    U = decode_type(jsondict[type_key])
    U <: type || U >: type || error("Type in file doesn't match target type: $(dict[type_key]) not a subtype of $T")

    if serialize_with_params(type)
      if isnothing(params) 
        params = load_type_params(state, type, jsondict[type_key][:params])
      end
      loaded = load_object(state, type, jsondict[:data], params)
    else
      Base.issingletontype(type) && return type()
      loaded = load_object(state, type, jsondict[:data])
    end
  else
    loaded = load_typed_object(state, jsondict; override_params=params)
  end

  if haskey(jsondict, :id)
    global_serializer_state.obj_to_id[loaded] = UUID(jsondict[:id])
    global_serializer_state.id_to_obj[UUID(jsondict[:id])] = loaded
  end

  return loaded
end

function load(filename::String; params::Any = nothing, type::Any = nothing)
  open(filename) do file
    return load(file; params=params, type=type)
  end
end
