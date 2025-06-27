using JSON3

################################################################################
# Type Serializers (converting types to strings)
convert_type_to_string(T::DataType) = sprint(show, T; context=:module=>Oscar)

################################################################################
# Serializers
abstract type OscarSerializer end

struct JSONSerializer <: OscarSerializer end
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
end

function begin_node(s::SerializerState)
  if !s.new_level_entry
    write(s.io, ",")
  else
    s.new_level_entry = false
  end
  if !isnothing(s.key)
    key = string(s.key)
    write(s.io, "\"$key\":")
    s.key = nothing
  end
end

function begin_dict_node(s::SerializerState)
  begin_node(s)
  write(s.io, "{")
end

function end_dict_node(s::SerializerState)
  write(s.io, "}")

  if s.new_level_entry
    # makes sure that entries after empty dicts add comma
    s.new_level_entry = false
  end
end

function begin_array_node(s::SerializerState)
  begin_node(s)
  write(s.io, "[")
end

function end_array_node(s::SerializerState)
  write(s.io, "]")

  if s.new_level_entry
    # makes sure that entries after empty arrays add comma
    s.new_level_entry = false
  end
end

function serialize_dict(f::Function, s::SerializerState)
  begin_dict_node(s)
  f()
  end_dict_node(s)
end

function serialize_array(f::Function, s::SerializerState)
  begin_array_node(s)
  f()
  end_array_node(s)
end

function set_key(s::SerializerState, key::Symbol)
  @req isnothing(s.key) "Key :$(s.key) is being overridden by :$key before write."
  s.key = key
end

## operations for an in-order tree traversals
## all nodes (dicts or arrays) contain all child nodes

function save_data_dict(f::Function, s::SerializerState,
                        key::Union{Symbol, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  serialize_dict(s) do
    s.new_level_entry = true
    f()
  end
end

function save_data_array(f::Function, s::SerializerState,
                         key::Union{Symbol, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  serialize_array(s) do
    s.new_level_entry = true
    f()
  end
end

function save_data_basic(s::SerializerState, x::Any,
                         key::Union{Symbol, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  begin_node(s)
  str = string(x)
  JSON.show_string(s.io, str)
  nothing
end

function save_data_json(s::SerializerState, jsonstr::Any,
                        key::Union{Symbol, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  begin_node(s)
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
  if !isempty(s.refs) 
    save_data_dict(s, refs_key) do
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

function handle_refs(s::SerializerState{IPCSerializer}) end

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
  obj::Union{Dict{Symbol, Any}, Vector, JSON3.Object, JSON3.Array, BasicTypeUnion}
  key::Union{Symbol, Int, Nothing}
  refs::Union{Dict{Symbol, Any}, JSON3.Object, Nothing}
  with_attrs::Bool
end

# general loading of a reference

function load_ref(s::DeserializerState)
  id = s.obj
  if haskey(global_serializer_state.id_to_obj, UUID(id))
    loaded_ref = global_serializer_state.id_to_obj[UUID(id)]
  else
    s.obj = s.refs[Symbol(id)]
    loaded_ref = load_typed_object(s)
    global_serializer_state.id_to_obj[UUID(id)] = loaded_ref
    global_serializer_state.obj_to_id[loaded_ref] = UUID(id)
  end
  return loaded_ref
end

function Base.haskey(s::DeserializerState, key::Symbol)
  s.obj isa String && return false
  load_node(s) do obj
    key in keys(obj)
  end
end

function set_key(s::DeserializerState, key::Union{Symbol, Int})
  @req isnothing(s.key) "Object at Key :$(s.key) hasn't been deserialized yet."
  s.key = key
end

function load_node(f::Function, s::DeserializerState,
                   key::Union{Symbol, Int, Nothing} = nothing)
  if s.obj isa String && !isnothing(tryparse(UUID, s.obj))
    return load_ref(s)
  end

  !isnothing(key) && set_key(s, key)
  obj = s.obj
  s.obj = isnothing(s.key) ? s.obj : s.obj[s.key]
  s.key = nothing
  result = f(s.obj)
  s.obj = obj
  return result
end

function load_array_node(f::Function, s::DeserializerState,
                         key::Union{Symbol, Int, Nothing} = nothing)
  load_node(s, key) do array
    [load_node(x -> f((i, x)), s, i) for (i, _) in enumerate(array)]
  end
end

function serializer_open(
  io::IO,
  serializer::OscarSerializer,
  with_attrs::Bool)
  
  # some level of handling should be done here at a later date
  return SerializerState(serializer, true, UUID[], io, nothing, with_attrs)
end

function deserializer_open(io::IO, serializer::OscarSerializer, with_attrs::Bool)
  obj = JSON3.read(io)
  refs = nothing
  if haskey(obj, refs_key)
    refs = obj[refs_key]
  end
  
  return DeserializerState(serializer, obj, nothing, refs, with_attrs)
end

function deserializer_open(io::IO, serializer::IPCSerializer, with_attrs::Bool) 
  # Using a JSON3.Object from JSON3 version 1.13.2 causes
  # put_params to hang
  #obj = JSON3.read(io)
  obj = JSON.parse(io, dicttype=Dict{Symbol, Any})

  return DeserializerState(serializer, obj, nothing, nothing, with_attrs)
end
