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
mutable struct SerializerState
  # dict to track already serialized objects
  new_level_entry::Bool
  # UUIDs that point to the objs in the global state,
  # ideally this would be an ordered set
  refs::Vector{UUID} 
  io::IO
  key::Union{Symbol, Nothing}
end

function SerializerState(io::IO)
  return SerializerState(true, UUID[], io, nothing)
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

struct DeserializerState
  # or perhaps Dict{Int,Any} to be resilient against corrupts/malicious files using huge ids
  # the values of refs are objects to be deserialized
  refs::Dict{Symbol, Dict}
end

function DeserializerState()
  return DeserializerState(Dict{Symbol, Any}())
end

function finish_writing(s::SerializerState)
  # nothing to do here
end

function set_key(s::SerializerState, key::Symbol)
  @req isnothing(s.key) "Key :$(s.key) is being overriden by :$key before write."
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
  write(s.io, "\"$str\"")
end

function save_data_json(s::SerializerState, jsonstr::Any,
                        key::Union{Symbol, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  begin_node(s)
  write(s.io, jsonstr)
end

function serializer_open(io::IO)
  # some level of handling should be done here at a later date
  return SerializerState(io)
end

function serializer_close(s::SerializerState)
  finish_writing(s)
end

function deserializer_open(io::IO)
  # should eventually take io
  return DeserializerState()
end
