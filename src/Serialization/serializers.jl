## JSON

mutable struct JSONSerializer
    open_objects::Array{Any}
end

function JSONSerializer()
    return JSONSerializer(Any[])
end

function serialize_dict(f::Function, s::JSONSerializer, key::Union{Symbol,Nothing} = nothing)
    begin_dict_node(s)
    f()
    end_node(s, key)
end

function serialize_array(f::Function, s::JSONSerializer, key::Union{Symbol,Nothing} = nothing)
    begin_array_node(s)
    f()
    end_node(s, key)
end

function begin_dict_node(s::JSONSerializer)
    push!(s.open_objects, Dict{Symbol, Any}())
end

function begin_array_node(s::JSONSerializer)
    push!(s.open_objects, Any[])
end

function end_node(s::JSONSerializer, key::Union{Symbol,Nothing} = nothing)
    obj = pop!(s.open_objects)
    add_object(s, obj, key)
end

function add_object(s::JSONSerializer, obj::Any, key::Symbol)
    s.open_objects[end][key] = obj
end

function add_object(s::JSONSerializer, obj::Any, key::Nothing=nothing)
    push!(s.open_objects[end], obj)
end

function finish_writing(s::JSONSerializer)
    @req length(s.open_objects) == 1 "too many open objects in serializer"
    str = json(s.open_objects[1])
    write(s.io, str)
end

## streaming

mutable struct StreamSerializer
    
end

function serialize_dict(f::Function, s::StreamSerializer, key::Union{Symbol,Nothing} = nothing)
    begin_dict_node(s, key)
    f()
    end_node(s)
end

function finish_writing(s::StreamSerializer)
    # nothing to do here
end

################################################################################
# (de)Serializer States

# struct which tracks state for (de)serialization
mutable struct SerializerState
    # dict to track already serialized objects
    objmap::IdDict{Any, UUID}
    depth::Int
    refs::Dict{Symbol, Any}
    serializer::JSONSerializer 
    key::Union{Symbol, Nothing}
    # TODO: if we don't want to produce intermediate dictionaries (which is a lot of overhead), we would probably store an `IO` object here
    # io::IO
end

function SerializerState()
    return SerializerState(IdDict{Any, UUID}(), 0, Dict(), JSONSerializer(), nothing)
end

struct DeserializerState
    objs::Dict{UUID, Any}  # or perhaps Dict{Int,Any} to be resilient against corrupts/malicious files using huge ids
    refs::Dict{Symbol, Any}
end

function DeserializerState()
    return DeserializerState(Dict{UUID, Any}(), Dict{Symbol,Any}())
end

## operations for an in-order tree traversals
## all nodes (dicts or arrays) contain all child nodes

function data_dict(f::Function, s::SerializerState)
    key = s.key
    s.key = nothing
    serialize_dict(s.serializer, key) do
        s.depth += 1
        f()
        s.depth -= 1
    end
end

function data_array(f::Function, s::SerializerState)
    key = s.key
    s.key = nothing
    serialize_array(s.serializer, key) do
        s.depth += 1
        f()
        s.depth -= 1
    end
end
function data_basic(s::SerializerState, x::Any)
    key = s.key
    s.key = nothing
    add_object(s.serializer, x, key)
end

function serializer_open(io::IO)
    # some level of handling should be done here at a later date
    return SerializerState()
end

function serializer_close(s::SerializerState)
    finish_writing(s.serializer)
end

#{
#  "key": [ 
#     {...},
#     {...}
#  ]
#}
