# this file kept as a reference  not to be used in
register_serializer(::Type{SomePolyRing})
register_serializer(::Type{SomePolyRingElem})

# not really correct, we need this for the object
type_needs_parent(::Type{SomePolyRingElem}) = true

function save_object(s::SerializerState, x::SomePolyRing)
   data_dict(s) do
      save_typed_object(s, base_ring(s), :base_ring)
      save_object(s, symbols(s), :symbols)
   end
end

function save_object(s::SerializerState, x::SomePolyRingElem)
   save_array(s, terms(x))
end

function save_array(s::SerializerState, x::AbstractArray)
   data_array(s) do
      for elem in x
         save_object(s, elem)
      end
   end
end

function save_object(s::SerializerState, x::String)
   data_basic(s, x)
end

# {
#   "type": { 
#     "name": "MPolyRingElem",
#     "parents": [
#         "fe6c33d9-ce7f-4f28-b650-597a5c87c887",
#         "23f25330-83f7-43a0-ac74-da6f2caa7eb8",
#         "a7029443-b1d3-4708-a66d-f68eb6616fcf"
#     ],
#   },
#   "data": [[[4, 3], [[0, "2"]]],
#           [[2, 0], [[0, "3"], [1, "1"]]],
#           [[0, 1], [[1, "5"]]],
#           [[0, 0], [[0, "1"]]]]
# }

function save_object(s::SerializerState, x::Any, key::Symbol)
   s.key=key
   save_object(s, x)
end


function save_typed_object(s::SerializerState, x::Any)
   if type_needs_parents(x)
      save_ref_type(s, parent(x), :type)
      save_object(s, x, :data)
   elseif is_singleton_type(x)
      save_object(s, encode_type(x), :type)
   elseif x isa AbstractVector
      save_object(s, "Vector", :type)
      save_object(s, innermost_elem_type(x), :elem_type)
      save_array(s, x, :data)
   else
      save_object(s, encode_type(x), :type)
      save_object(s, x, :data)
   end
end


function save(io::IO, obj::Any, metadata::Any=nothing)
   state = serializer_open(io)
   data_dict(s) do
      add_header(s, metadata)
      save_typed_object(s, obj)
   end
   serializer_close(state)
end
