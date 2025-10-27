module Serialization

using ..Oscar
using UUIDs
import JSON

using ..Oscar: _grading,
  FreeAssociativeAlgebraIdeal,
  IdealGens,
  LaurentMPolyIdeal,
  MPolyAnyMap,
  MPolyLocalizedRingHom,
  MPolyQuoLocalizedRingHom,
  NormalToricVarietyType,
  Orderings,
  PhylogeneticTree,
  pm_object,
  PolyhedralObject,
  scalar_types,
  VERSION_NUMBER,
  _pmdata_for_oscar

using ..Oscar: is_terse, Lowercase, pretty, terse

using Distributed: RemoteChannel

# This type should not be exported and should be before serializers
const BasicTypeUnion = Union{String, QQFieldElem, Symbol,
                       Number, ZZRingElem, TropicalSemiringElem}

include("serializers.jl")

const type_key = :_type

################################################################################
# Meta Data

@Base.kwdef struct MetaData
  author_orcid::Union{String, Nothing} = nothing
  name::Union{String, Nothing} = nothing
  description::Union{String, Nothing} = nothing
end

# FIXME: this function is exported but undocumented
function metadata(;args...)
  return MetaData(;args...)
end

# FIXME: this function is exported but undocumented
function read_metadata(filename::String)
  open(filename) do io
    obj = JSON3.read(io)
    println(JSON.json(obj[:meta], 2))
  end
end

################################################################################
# Serialization info

function serialization_version_info(obj::AbstractDict{Symbol, Any})
  ns = obj[:_ns]
  version_info = ns[:Oscar][2]
  return version_number(version_info)
end

function version_number(v_number::String)
  return VersionNumber(v_number)
end

# needed for older versions
function version_number(dict::AbstractDict)
  return VersionNumber(dict[:major], dict[:minor], dict[:patch])
end

const oscar_serialization_version = Ref{Dict{Symbol, Any}}()

function get_oscar_serialization_version()
  if isassigned(oscar_serialization_version)
    return oscar_serialization_version[]
  end
  if Oscar.is_dev
    next_version = "$(VERSION_NUMBER.major).$(VERSION_NUMBER.minor).$(VERSION_NUMBER.patch)"
    n_upgrades = count(x -> startswith(x, next_version) && endswith(x, ".jl"),
                       readdir(joinpath(@__DIR__, "Upgrades")))


    commit_hash = get(Oscar._get_oscar_git_info(), :commit, "unknown")
    version_info = iszero(n_upgrades) ? "$VERSION_NUMBER-$commit_hash" : "$VERSION_NUMBER-$n_upgrades-$commit_hash"
    result = Dict{Symbol, Any}(
      :Oscar => ["https://github.com/oscar-system/Oscar.jl", version_info]
    )
  else
    result = Dict{Symbol, Any}(
      :Oscar => ["https://github.com/oscar-system/Oscar.jl", VERSION_NUMBER]
    )
  end
  return oscar_serialization_version[] = result
end

################################################################################
# Type attribute map
const type_attr_map = Dict{String, Vector{Symbol}}()

attrs_list(T::Type) = get(type_attr_map, encode_type(T), Symbol[])

with_attrs(s::T) where T <: Union{DeserializerState, SerializerState} = s.with_attrs

################################################################################
# (De|En)coding types

# parameters of type should not matter here
const reverse_type_map = Dict{String, Union{Dict{String, Type}, Type}}()

function encode_type(::Type{T}) where T
  error(
    """Unsupported type '$T' for encoding. To add support see
    https://docs.oscar-system.org/stable/DeveloperDocumentation/serialization/
    """
  )
end

function decode_type(s::String)
  return get(reverse_type_map, s) do
    error("unsupported type '$s' for decoding")
  end
end

function decode_type(s::DeserializerState)
  if s.obj isa String
    if !isnothing(tryparse(UUID, s.obj))
      id = s.obj
      obj = s.obj
      if isnothing(s.refs)
        return typeof(global_serializer_state.id_to_obj[UUID(id)])
      end
      s.obj = s.refs[Symbol(id)]
      T = decode_type(s)
      s.obj = obj
      return T
    end
    return decode_type(s.obj)
  end

  if type_key in keys(s.obj)
    return load_node(s, type_key) do _
      decode_type(s)
    end
  end

  if :name in keys(s.obj)
    if :_instance in keys(s.obj)
      return get(reverse_type_map[s.obj[:name]], s.obj[:_instance]) do
        unsupported_instance = s.obj[:_instance]
        error("unsupported instance '$unsupported_instance' for decoding")
      end
    else
      return load_node(s, :name) do _
        decode_type(s)
      end
    end
  end
end

################################################################################
# TypeParams Struct
struct TypeParams{T, S}
  type::Type{T}
  params::S

  function TypeParams(T::Type, args::Pair...)
    return new{T, typeof(args)}(T, args)
  end
  TypeParams(T::Type, obj) = new{T, typeof(obj)}(T, obj)
end

params(tp::TypeParams) = tp.params
type(tp::TypeParams) = tp.type

type_params(obj::T) where T = TypeParams(T, nothing)

function Base.show(io::IO, tp::TypeParams{T, Tuple}) where T
  if is_terse(io)
    print(io, "Type parameters for $T")
  else
    io = pretty(io)
    print(io, "Type parameters for $T")
    for param in params(tp)
      println(io, "")
      print(terse(io), Lowercase(), param)
    end
  end
end

function Base.show(io::IO, tp::TypeParams{T, S}) where {T, S}
  if is_terse(io)
    print(io, "Type parameters for $T")
  else
    io = pretty(io)
    print(io, "Type parameters for $T ")
    print(terse(io), Lowercase(), params(tp))
  end
end

# ATTENTION
# We need to distinguish between data with a globally defined normal form and data where such a normal form depends on some parameters.
# In particular, this does NOT ONLY depend on the type; see, e.g., FqField.

################################################################################
# High level

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
  if Base.issingletontype(T)
    save_object(s, encode_type(T), type_key)
  else
    save_type_params(s, x, type_key)
    save_object(s, x, :data)
  end
  if with_attrs(s)
    attrs = attrs_list(T)
    !isempty(attrs) && save_attrs(s, x)
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

################################################################################
# (save | load) TypeParams

function save_type_params(s::SerializerState, obj::Any, key::Symbol)
  set_key(s, key)
  save_type_params(s, obj)
end

function save_type_params(s::SerializerState, obj::T) where T
  save_type_params(s, type_params(obj))
end

function save_type_params(s::SerializerState, tp::TypeParams)
  save_data_dict(s) do
    T = type(tp)
    type_encoding = encode_type(T)
    if reverse_type_map[type_encoding] isa Dict
      save_object(s, convert_type_to_string(T), :_instance)
    end
    
    save_object(s, type_encoding, :name)
    # params(tp) isa TypeParams if the type isa container type
    if params(tp) isa TypeParams
      save_type_params(s, params(tp), :params)
    else
      save_typed_object(s, params(tp), :params)
    end
  end
end

function save_type_params(s::SerializerState,
                          ::TypeParams{T, Nothing}) where T
  type_encoding = encode_type(T)
  if reverse_type_map[type_encoding] isa Dict
    save_data_dict(s) do
      save_object(s, type_encoding, :name)
      save_object(s, convert_type_to_string(T), :_instance)
    end
  else
    save_object(s, type_encoding)
  end
end

function save_type_params(s::SerializerState,
                          tp::TypeParams{<:TypeParams, <:Tuple{Vararg{Pair}}})
  for param in params(tp)
    save_type_params(s, param.second, Symbol(param.first))
  end
end

function save_type_params(s::SerializerState,
                          tp::TypeParams{<:TypeParams, <:Tuple})
  save_data_array(s) do 
    for param in params(tp)
      save_type_params(s, param)
    end
  end
end

function save_type_params(s::SerializerState,
                          tp::TypeParams{T, <:Tuple{Vararg{Pair}}}) where T
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_dict(s, :params) do
      for param in params(tp)
        if param.second isa Type
          save_object(s, encode_type(param.second), Symbol(param.first))
        elseif !(param.second isa TypeParams)
          if param.second isa Tuple
            save_data_array(s, Symbol(param.first)) do
              for entry in param.second
                if serialize_with_id(entry)
                  save_object(s, save_as_ref(s, entry))
                else
                  save_data_dict(s) do
                    save_typed_object(s, entry)
                  end
                end
              end
            end
          else
            save_typed_object(s, param.second, Symbol(param.first))
          end
        else
          save_type_params(s, param.second, Symbol(param.first))
        end
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type, key::Symbol)
  load_node(s, key) do _
    load_type_params(s, T)
  end
end

function load_type_array_params(s::DeserializerState)
  load_array_node(s) do obj
    T = decode_type(s)
    if obj isa String
      !isnothing(tryparse(UUID, s.obj)) && return load_ref(s)
      return T
    end
    return load_type_params(s, T)[2]
  end
end

function load_type_params(s::DeserializerState, T::Type)
  if s.obj isa String
    if !isnothing(tryparse(UUID, s.obj))
      return T, load_ref(s)
    end
    return T, nothing
  end
  if haskey(s, :params)
    load_node(s, :params) do obj
      if obj isa JSON3.Array || obj isa Vector
        params = load_type_array_params(s)
      elseif obj isa String || haskey(s, :params)
        U = decode_type(s)
        if Base.issingletontype(U)
          params = U()
        else
          params = load_type_params(s, U)[2]
        end
      # handle cases where type_params is a dict of params
      elseif !haskey(obj, type_key) 
        params = Dict{Symbol, Any}()
        for (k, _) in obj
          params[k] = load_node(s, k) do obj
            if obj isa JSON3.Array || obj isa Vector
              return load_type_array_params(s)
            end
            
            U = decode_type(s)
            if obj isa String && isnothing(tryparse(UUID, obj))
              return U
            end
            return load_type_params(s, U)[2]
          end
        end
      else
        params = load_typed_object(s)
      end
      # all types where the type T should be updated with a subtype i.e. T -> T{U}
      # need to implement their own method, see for example containers
      return T, params
    end
  elseif haskey(s, :_instance)
    T, nothing
  else
    return T, load_typed_object(s)
  end
end

function load_typed_object(s::DeserializerState, key::Symbol; override_params::Any = nothing)
  load_node(s, key) do _
    load_typed_object(s; override_params=override_params)
  end
end

# The load mechanism first checks if the type needs to load necessary
# parameters before loading it's data, if so a type tree is traversed
function load_typed_object(s::DeserializerState; override_params::Any = nothing)
  T = decode_type(s)
  Base.issingletontype(T) && return T()
  if !isnothing(override_params)
    T, _ = load_type_params(s, T, type_key)
    params = override_params
  else
    s.obj isa String && !isnothing(tryparse(UUID, s.obj)) && return load_ref(s)
    T, params = load_type_params(s, T, type_key)
  end
  obj = load_node(s, :data) do _
    return load_object(s, T, params)
  end
  load_attrs(s, obj)
  return obj
end

function load_object(s::DeserializerState, T::Type, key::Union{Symbol, Int})
  load_node(s, key) do _
    load_object(s, T)
  end
end

function load_object(s::DeserializerState, T::Type, params::S,
                     key::Union{Symbol, Int}) where S
  load_node(s, key) do _
    load_object(s, T, params)
  end
end

load_object(s::DeserializerState, T::Type, ::Nothing) = load_object(s, T)

################################################################################
# serializing attributes
function save_attrs(s::SerializerState, obj::T) where T
  if any(attr -> has_attribute(obj, attr), attrs_list(T))
    save_data_dict(s, :attrs) do
      for attr in attrs_list(T)
        has_attribute(obj, attr) && save_typed_object(s, get_attribute(obj, attr), attr)
      end
    end
  end
end

function load_attrs(s::DeserializerState, obj::T) where T
  !with_attrs(s) && return

  haskey(s, :attrs) && load_node(s, :attrs) do d
    for attr in keys(d)
      set_attribute!(obj, attr, load_typed_object(s, attr))
    end
  end
end

################################################################################
# Type Registration
function register_serialization_type(@nospecialize(T::Type), str::String)
  if haskey(reverse_type_map, str) 
    init = reverse_type_map[str]
    # promote the value to a dictionary if necessary
    if init isa Type
      init = Dict{String, Type}(convert_type_to_string(init) => init)
    end
    reverse_type_map[str] = merge(Dict{String, Type}(convert_type_to_string(T) => T), init)
  else
    reverse_type_map[str] = T
  end
end

function register_attr_list(@nospecialize(T::Type),
                            attrs::Union{Vector{Symbol}, Nothing})
  if !isnothing(attrs)
    serialize_with_id(T) || error("Only types that are stored as references can store attributes")
    type_attr_map[encode_type(T)] = attrs
  end
end

import Serialization.serialize
import Serialization.deserialize
import Serialization.serialize_type
import Distributed.AbstractSerializer

# add these here so that the proper errors are thrown
# when the type hasn't been registered
serialize_with_id(::Type) = false
serialize_with_id(obj::Any) = false

function register_serialization_type(ex::Any, str::String, uses_id::Bool, attrs::Any)
  return esc(
    quote
      Oscar.Serialization.register_serialization_type($ex, $str)
      Oscar.Serialization.encode_type(::Type{<:$ex}) = $str
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

      Oscar.Serialization.serialize_with_id(obj::T) where T <: $ex = $uses_id
      Oscar.Serialization.serialize_with_id(T::Type{<:$ex}) = $uses_id

      # add list of possible attributes to save for a given type to a global dict
      Oscar.Serialization.register_attr_list($ex, $attrs)

      # only extend serialize on non std julia types
      if !($ex <: Union{Number, String, Bool, Symbol, Vector, Tuple, Matrix, NamedTuple, Dict, Set, Array})
        function Oscar.Serialization.serialize(s::Oscar.Serialization.AbstractSerializer, obj::T) where T <: $ex
          Oscar.Serialization.serialize_type(s, T)
          Oscar.Serialization.save(s.io, obj; serializer=Oscar.Serialization.IPCSerializer())

        end
        function Oscar.Serialization.deserialize(s::Oscar.Serialization.AbstractSerializer, T::Type{<:$ex})
          Oscar.Serialization.load(s.io; serializer=Oscar.Serialization.IPCSerializer())
        end
      end
    end)
end

"""
    @register_serialization_type NewType "String Representation of type" uses_id uses_params [:attr1, :attr2]

`@register_serialization_type` is a macro to ensure that the string we generate
matches exactly the expression passed as first argument, and does not change
in unexpected ways when import/export statements are adjusted.

Passing a string argument will override how the type is stored as a string.

When setting `uses_id` the object will be stored as a reference and
will be referred to throughout the serialization sessions using a `UUID`.
This should typically only be used for types that do not have a fixed
normal form for example `PolyRing` and `MPolyRing`.

Using the `uses_params` flag will serialize the object with a more structured type
description which will make the serialization more efficient see the discussion on
`save_type_params` / `load_type_params` below.

Passing a vector of symbols that correspond to attributes of type
indicates which attributes will be serialized when using save with `with_attrs=true`.

"""
macro register_serialization_type(ex::Any, args...)
  uses_id = false
  str = nothing
  attrs = nothing
  for el in args
    if el isa String
      str = el
    elseif el == :uses_id
      uses_id = true
    else
      attrs = el
    end
  end
  if str === nothing
    # here we use string since on an expression.
    # this choice means we should write types without the namespace in front
    # when registering, and in convert_type_to_string
    str = string(ex)
  end

  return register_serialization_type(ex, str, uses_id, attrs)
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
include("MPolyMap.jl")
include("Algebras.jl")
include("polymake.jl")
include("TropicalGeometry.jl")
include("QuadForm.jl")
include("GAP.jl")
include("Groups.jl")
include("LieTheory.jl")

include("Upgrades/main.jl")
include("parallel.jl")


################################################################################
# Interacting with IO streams and files

"""
    save(io::IO, obj::Any; metadata::MetaData=nothing, with_attrs::Bool=true)
    save(filename::String, obj::Any, metadata::MetaData=nothing, with_attrs::Bool=true)

Save an object `obj` to the given io stream
respectively to the file `filename`. When used with `with_attrs=true` then the object will
save it's attributes along with all the attributes of the types used in the object's struct.
The attributes that will be saved are defined during type registration see
[`@register_serialization_type`](@ref)

See [`load`](@ref).

# Examples

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current))
julia> meta = metadata(author_orcid="0000-0000-0000-0042", name="42", description="The meaning of life, the universe and everything")
Oscar.Serialization.MetaData("0000-0000-0000-0042", "42", "The meaning of life, the universe and everything")

julia> save("fourtitwo.mrdi", 42; metadata=meta);

julia> read_metadata("fourtitwo.mrdi")
{
  "author_orcid": "0000-0000-0000-0042",
  "name": "42",
  "description": "The meaning of life, the universe and everything"
}

julia> load("fourtitwo.mrdi")
42
```
"""
function save(io::IO, obj::T; metadata::Union{MetaData, Nothing}=nothing,
              with_attrs::Bool=true,
              serializer::OscarSerializer = JSONSerializer()) where T
  s = serializer_open(io, serializer, with_attrs)
  save_data_dict(s) do 
    # write out the namespace first
    save_header(s, get_oscar_serialization_version(), :_ns)

    save_typed_object(s, obj)

    if serialize_with_id(T)
      ref = get(global_serializer_state.obj_to_id, obj, nothing)
      if isnothing(ref)
        ref = global_serializer_state.obj_to_id[obj] = uuid4()
        global_serializer_state.id_to_obj[ref] = obj
      end
      save_object(s, string(ref), :id)
    end

    handle_refs(s)

    if !isnothing(metadata)
      save_json(s, JSON3.write(metadata), :meta)
    end
  end
  serializer_close(s)
  return nothing
end

function save(filename::String, obj::Any;
              metadata::Union{MetaData, Nothing}=nothing,
              serializer::OscarSerializer=JSONSerializer(),
              with_attrs::Bool=true)
  dir_name = dirname(filename)
  # julia dirname does not return "." for plain filenames without any slashes
  temp_file = tempname(isempty(dir_name) ? pwd() : dir_name)
  
  open(temp_file, "w") do file
    save(file, obj;
         metadata=metadata,
         with_attrs=with_attrs,
         serializer=serializer)
  end
  Base.Filesystem.rename(temp_file, filename) # atomic "multi process safe"
  return nothing
end

"""
    load(io::IO; params::Any = nothing, type::Any = nothing, with_attrs::Bool=true)
    load(filename::String; params::Any = nothing, type::Any = nothing, with_attrs::Bool=true)

Load the object stored in the given io stream
respectively in the file `filename`.

If `params` is specified, then the root object of the loaded data
either will attempt a load using these parameters. In the case of Rings this
results in setting its parent, or in the case of a container of ring types such as
`Vector` or `Tuple`, then the parent of the entries will be set using their
 `params`.

If a type `T` is given then attempt to load the root object of the data
being loaded with this type; if this fails, an error is thrown.

If `with_attrs=true` the object will be loaded with attributes available from
the file (or serialized data).

See [`save`](@ref).

# Examples

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current))
julia> save("fourtitwo.mrdi", 42);

julia> load("fourtitwo.mrdi")
42

julia> load("fourtitwo.mrdi"; type=Int64)
42

julia> R, x = QQ[:x]
(Univariate polynomial ring in x over QQ, x)

julia> p = x^2 - x + 1
x^2 - x + 1

julia> save("p.mrdi", p)

julia> p_loaded = load("p.mrdi", params=R)
x^2 - x + 1

julia> parent(p_loaded) === R
true

julia> save("p_v.mrdi", [p, p])

julia> loaded_p_v = load("p_v.mrdi", params=R)
2-element Vector{QQPolyRingElem}:
 x^2 - x + 1
 x^2 - x + 1

julia> parent(loaded_p_v[1]) === parent(loaded_p_v[2]) === R
true
```
"""
function load(io::IO; params::Any = nothing, type::Any = nothing,
              serializer=JSONSerializer(), with_attrs::Bool=true)
  s = deserializer_open(io, serializer, with_attrs)
  if haskey(s.obj, :id)
    id = s.obj[:id]
    if haskey(global_serializer_state.id_to_obj, UUID(id))
      return global_serializer_state.id_to_obj[UUID(id)]
    end
  end

  # handle different namespaces
  polymake_obj = load_node(s) do d
    @req :_ns in keys(d) "Namespace is missing"
    load_node(s, :_ns) do _ns
      if :polymake in keys(_ns)
        return load_from_polymake(Dict(d))
      end
    end
  end
  if !isnothing(polymake_obj)
    return polymake_obj
  end

  load_node(s, :_ns) do _ns
    @req haskey(_ns, :Oscar) "Not an Oscar object"
  end

  # deal with upgrades
  file_version = load_node(s) do obj
    serialization_version_info(obj)
  end
  if file_version < VERSION_NUMBER
    # we need a mutable dictionary
    jsondict = copy(s.obj)
    jsondict = upgrade(file_version, jsondict)
    jsondict_str = JSON3.write(jsondict)
    s = deserializer_open(IOBuffer(jsondict_str),
                          serializer,
                          with_attrs)
  end
  
  try
    if params isa TypeParams
      params = _convert_override_params(params)
    end
    if type !== nothing
      # Decode the stored type, and compare it to the type `T` supplied by the caller.
      # If they are identical, just proceed. If not, then we assume that either
      # `T` is concrete, in which case `T <: U` should hold; or else `U` is
      # concrete, and `U <: T` should hold.
      #
      # This check should maybe change to a check on the whole type tree?
      U = load_node(s, type_key) do _
        decode_type(s)
      end
      
      U <: type || U >: type || error("Type in file doesn't match target type: $(dict[type_key]) not a subtype of $type")

      Base.issingletontype(type) && return type()
      if isnothing(params)
        _, params = load_node(s, type_key) do _
          load_type_params(s, U)
        end
      end
      load_node(s, :data) do _
        loaded = load_object(s, type, params)
      end
    else
      loaded = load_typed_object(s; override_params=params)
    end

    if :id in keys(s.obj)
      load_node(s, :id) do id
        global_serializer_state.obj_to_id[loaded] = UUID(id)
        global_serializer_state.id_to_obj[UUID(id)] = loaded
      end
    end
    return loaded
  catch e
    if VersionNumber(replace(string(file_version), r"DEV.+" => "DEV")) > VERSION_NUMBER
      @warn """
      Attempted loading file stored with Oscar version $file_version
      using Oscar version $VERSION_NUMBER
      """
    end

    if contains(string(file_version), "DEV")
      commit = split(string(file_version), "-")[end]
      @warn "Attempted loading file stored using a DEV version with commit $commit"
    end
    rethrow(e)
  end
end

function load(filename::String; params::Any = nothing,
              type::Any = nothing, with_attrs::Bool=true,
              serializer::OscarSerializer=JSONSerializer())
  open(filename) do file
    return load(file; params=params, type=type, serializer=serializer)
  end
end

_convert_override_params(tp::TypeParams{T, S}) where {T, S} = _convert_override_params(params(tp))
_convert_override_params(tp::TypeParams{T, <:Tuple{Vararg{Pair}}}) where T = Dict(_convert_override_params(params(tp)))
_convert_override_params(tp::TypeParams{T, S}) where {T <: MatVecType, S} = _convert_override_params(params(tp))
_convert_override_params(tp::TypeParams{T, S}) where {T <: Set, S} = _convert_override_params(params(tp))
_convert_override_params(tp::TypeParams{<: NamedTuple, S}) where S = _convert_override_params(values(params(tp)))
_convert_override_params(tp::TypeParams{T, Nothing}) where T <: Oscar.AbstractGraph = T
_convert_override_params(tp::TypeParams{<:Array, <:Tuple{Vararg{Pair}}}) = Dict(_convert_override_params(params(tp)))[:subtype_params]

function _convert_override_params(tp::TypeParams{<:Dict, <:Tuple{Vararg{Pair}}})
  vp_pair = filter(x -> :value_params == x.first, params(tp))
  kp_pair = filter(x -> :key_params == x.first, params(tp))
  if !isempty(vp_pair)
    ov_params = Dict(k => _convert_override_params(v) for (k, v) in params(tp))
    if type(first(kp_pair).second) <: Union{Symbol, String, Int}
      return _convert_override_params(first(vp_pair).second)
    end
  else
    return Dict(k => (type(v), _convert_override_params(v)) for (k, v) in params(tp))
  end
end

_convert_override_params(obj::Any) = obj

#handles empty tuple ambiguity
_convert_override_params(obj::Tuple{}) = ()

_convert_override_params(t::Tuple{Vararg{TypeParams}}) = map(_convert_override_params, t)

function _convert_override_params(t::Tuple{Vararg{Pair}})
  map(x -> x.first => _convert_override_params(x.second), t)
end

# handle special polyhedral case
function _convert_override_params(tp::TypeParams{<:PolyhedralObject, <:Tuple{Vararg{Pair}}})
  # special treatement for the polymake parameters
  poly_params = Dict()
  for (k, v) in params(tp)
    if k == :pm_params
      poly_params[k] = Dict()
      for (pm_k, pm_v) in params(v)
        poly_params[k][pm_k] = (type(pm_v), _convert_override_params(params(pm_v)))
      end
    else
      poly_params[k] = v
    end
  end
  return poly_params
end

# handle monomial ordering
_convert_override_params(tp::TypeParams{T, S}) where {T <: MonomialOrdering, S} = T


export @register_serialization_type
export DeserializerState
export encode_type
export load
export load_array_node
export load_attrs
export load_node
export load_object
export load_ref
export save
export save_as_ref
export save_attrs
export save_data_array
export save_data_basic
export save_data_dict
export save_data_json
export save_object
export SerializerState
export serialize_with_id
export set_key
export TypeParams
export type_params
export with_attrs

end # module Serialization

using Oscar.Serialization
import Oscar.Serialization: load_object, save_object, type_params
import Oscar.Serialization: reset_global_serializer_state

# FIXME: the following functions are exported by us but undocumented
import Oscar.Serialization: metadata, read_metadata
