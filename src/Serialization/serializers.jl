using JSON3

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
  write(s.io, "\"$str\"")
end

function save_data_json(s::SerializerState, jsonstr::Any,
                        key::Union{Symbol, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  begin_node(s)
  write(s.io, jsonstr)
end

function serializer_close(s::SerializerState)
  finish_writing(s)
end
  
function finish_writing(s::SerializerState)
  # nothing to do here
end

mutable struct DeserializerState
  # or perhaps Dict{Int,Any} to be resilient against corrupts/malicious files using huge ids
  # the values of refs are objects to be deserialized
  obj::Union{JSON3.Object, JSON3.Array, BasicTypeUnion}
  key::Union{Symbol, Nothing}
  refs::Union{JSON3.Object, Nothing}
end

function set_key(s::DeserializerState, key::Union{Symbol, Int})
  @req isnothing(s.key) "Object at Key :$(s.key) hasn't been deserialized yet."

  s.key = key
end


function deserialize_node(f::Function, s::DeserializerState)
  f(s.obj)
end

function load_node(f::Function, s::DeserializerState,
                    key::Union{Symbol, Int, Nothing} = nothing)
  !isnothing(key) && set_key(s, key)
  obj = deepcopy(s.obj)
  s.obj = isnothing(s.key) ? s.obj : s.obj[s.key]
  s.key = nothing
  result = f()
  s.obj = obj
  return result
end

################################################################################
# Serializers
abstract type OscarSerializer end

struct JSONSerializer <: OscarSerializer
  state::S where S <: Union{SerializerState, DeserializerState}
end

struct IPCSerializer <: OscarSerializer
  state::S where S <: Union{SerializerState, DeserializerState}
end

state(s::OscarSerializer) = s.state

function serializer_open(io::IO, T::Type{<: OscarSerializer})
  # some level of handling should be done here at a later date
  return T(SerializerState(true, UUID[], io, nothing))
end

function deserializer_open(io::IO, T::Type{<: OscarSerializer})
  obj = JSON3.read(io)
  refs = nothing
  if refs_key in keys(obj)
    refs = obj[refs_key]
  end
  
  return T(DeserializerState(obj, nothing, refs))
end
