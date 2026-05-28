################################################################################
# Type Serializers (converting types to strings)
convert_type_to_string(T::DataType) = sprint(show, T; context=:module=>Oscar)

################################################################################
# Serializers
abstract type OscarSerializer end

struct JSONSerializer <: OscarSerializer
  serialize_refs::Bool

  function JSONSerializer(; serialize_refs::Bool = true)
    return new(serialize_refs)
  end
end
struct IPCSerializer <: OscarSerializer end

abstract type MultiFileSerializer <: OscarSerializer end

struct LPSerializer <: MultiFileSerializer
  basepath::String
end

basepath(serializer::MultiFileSerializer) = serializer.basepath

################################################################################
# (de)Serializer States

# this struct is used to keep a global state of serialized objs during a session
# and allows sessions to store related objects across files
mutable struct GlobalSerializerState
  obj_to_id::IdDict{Any, UUID}
  id_to_obj::Dict{UUID, Any}
end

function GlobalSerializerState()
  return GlobalSerializerState(IdDict{Any, UUID}(), Dict{UUID, Any}())
end

const global_serializer_state = GlobalSerializerState()

function reset_global_serializer_state()
  empty!(global_serializer_state.obj_to_id)
  empty!(global_serializer_state.id_to_obj)
end

# struct which tracks state for (de)serialization
mutable struct SerializerState{T <: OscarSerializer}
  serializer::T
  new_level_entry::Bool
  # UUIDs that point to the objs in the global state,
  # ideally this would be an ordered set
  refs::Vector{UUID}
  io::IO
  key::Union{Symbol, Nothing}
  with_attrs::Bool
  pretty_print::Bool
end

function begin_node(s::SerializerState, key::Union{Symbol, Nothing})
  !isnothing(key) && set_key(s, key)
  if !s.new_level_entry
    if s.pretty_print
      println(s.io, ",")
    else
      write(s.io, ",")
    end
  else
    s.new_level_entry = false
  end
  if !isnothing(s.key)
    key = string(s.key)
    write(s.io, "\"$key\":")
    s.key = nothing
  end
end

function set_key(s::SerializerState, key::Symbol)
  @req isnothing(s.key) "Key :$(s.key) is being overridden by :$key before write."
  s.key = key
end

## operations for an in-order tree traversals
## all nodes (dicts or arrays) contain all child nodes

function _save_data_container(f::Function, s::SerializerState,
                        key::Union{Symbol, Nothing}, start::String, stop::String)
  begin_node(s, key)
  if s.pretty_print
    println(s.io, start)
    print(s.io, Indent())
  else
    write(s.io, start)
  end
  s.new_level_entry = true
  f()
  if s.pretty_print
    println(s.io, "")
    print(s.io, Dedent(), stop)
  else
    write(s.io, stop)
  end

  if s.new_level_entry
    # makes sure that entries after empty arrays or dicts add comma
    s.new_level_entry = false
  end
end

function save_data_dict(f::Function, s::SerializerState,
                        key::Union{Symbol, Nothing} = nothing)
  _save_data_container(f, s, key, "{", "}")
end

function save_data_array(f::Function, s::SerializerState,
                         key::Union{Symbol, Nothing} = nothing)
  _save_data_container(f, s, key, "[", "]")
end

function save_data_basic(s::SerializerState, x::Any,
                         key::Union{Symbol, Nothing} = nothing)
  begin_node(s, key)
  data = x isa Bool ? x : string(x)
  if s.pretty_print
    print(s.io, "")
    JSON.json(s.io, data)
    print(s.io, "")
  else
    JSON.json(s.io, data)
  end
  nothing
end

function save_data_json(s::SerializerState, jsonstr::Any,
                        key::Union{Symbol, Nothing} = nothing)
  begin_node(s, key)
  write(s.io, jsonstr)
end


function save_as_ref(s::SerializerState, obj::T) where T
  # find ref or create one
  ref = get(global_serializer_state.obj_to_id, obj, nothing)
  if !isnothing(ref)
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

function handle_refs(s::SerializerState)
  should_handle_refs(s.serializer) || return nothing
  if !isempty(s.refs) 
    save_data_dict(s, :_refs) do
      for id in s.refs
        ref_obj = global_serializer_state.id_to_obj[id]
        s.key = Symbol(id)
        save_data_dict(s) do
          save_typed_object(s, ref_obj)
        end
      end
    end
  end
end

should_handle_refs(::OscarSerializer) = true
should_handle_refs(s::JSONSerializer) = s.serialize_refs
should_handle_refs(::IPCSerializer) = false

function serializer_close(s::SerializerState)
  finish_writing(s)
end

function finish_writing(s::SerializerState)
  # nothing to do here
end

mutable struct DeserializerState{T <: OscarSerializer}
  # or perhaps Dict{Int,Any} to be resilient against corrupts/malicious files using huge ids
  # the values of refs are objects to be deserialized
  serializer::T
  obj::Union{JSON.LazyValue, BasicTypeUnion, Nothing}
  key::Union{Symbol, Int, Nothing}
  refs::Dict{UUID, JSON.LazyValues}
  with_attrs::Bool
end

function load_json(s::DeserializerState, ::Type{T}) where T
  return JSON.parse(s.obj, T)
end


# general loading of a reference

function load_ref(s::DeserializerState)
  id = load_json(s, UUID)
  if haskey(global_serializer_state.id_to_obj, id)
    loaded_ref = global_serializer_state.id_to_obj[id]
  else
    s.obj = s.refs[id]
    loaded_ref = load_typed_object(s)
    global_serializer_state.id_to_obj[id] = loaded_ref
    global_serializer_state.obj_to_id[loaded_ref] = id
  end
  return loaded_ref
end

function Base.isempty(s::DeserializerState)
  return iszero(length(s.obj))
end

function Base.haskey(s::DeserializerState, key::Symbol)
  !(s.obj isa JSON.LazyValue) && return false
  obj = s.obj[]
  obj isa String && return false
  return haskey(obj, key)
end

function set_key(s::DeserializerState, key::Union{Symbol, Int})
  @req isnothing(s.key) "Object at Key :$(s.key) hasn't been deserialized yet."
  s.key = key
end

function load_node(f::Function, s::DeserializerState,
                   key::Union{Symbol, Int, Nothing} = nothing)
  if is_string(s) && !isnothing(tryparse(UUID, load_json(s, String)))
    return load_ref(s)
  end

  !isnothing(key) && set_key(s, key)
  lazy_obj = s.obj
  if !isnothing(s.key)
    s.obj = s.obj[s.key]
  end
  s.key = nothing
  if isnothing(s.obj)
    result = nothing
  else
    result = f()
  end
  s.obj = lazy_obj
  return result
end

function load_array_node(f::Function, s::DeserializerState,
                         key::Union{Symbol, Int, Nothing} = nothing,
                         entry_type::DataType = Any)
  load_node(s, key) do
    i = 1
    result = entry_type[]
    sizehint(result, length(s.obj))
    foreach(s.obj) do v
      s.obj = v
      push!(result, f(i))
      i += 1
    end
    return result
  end
end

function serializer_open(
  io::IO,
  serializer::OscarSerializer,
  with_attrs::Bool,
  pretty_print::Bool)
  
  # some level of handling should be done here at a later date
  return SerializerState(serializer, true, UUID[], io, nothing, with_attrs, pretty_print)
end

function deserializer_open(io::IO, serializer::OscarSerializer, with_attrs::Bool)
  obj = JSON.lazy(io)
  refs_from_file = get(obj, :_refs, nothing)
  refs = Dict{UUID,JSON.LazyValues}()
  for k in propertynames(refs_from_file)
    refs[UUID(String(k))] = refs_from_file[k]
  end
  
  return DeserializerState(serializer, obj, nothing, refs, with_attrs)
end

function deserializer_open(io::IO, serializer::IPCSerializer, with_attrs::Bool)
  str = readuntil(io, '}'; keep=true)
  while !JSON.isvalidjson(str)
    str *= readuntil(io, '}'; keep=true)
  end
  obj = JSON.lazy(str)

  return DeserializerState(serializer, obj, nothing, Dict{UUID, JSON.LazyValues}(), with_attrs)
end
