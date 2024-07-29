@register_serialization_type AbstractLieAlgebra

function save_object(s::SerializerState, L::AbstractLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_object(s, _struct_consts(L), :struct_consts)
    save_object(s, symbols(L), :symbols)
  end
end

function load_object(s::DeserializerState, ::Type{AbstractLieAlgebra})
  R = load_typed_object(s, :base_ring)
  struct_consts = load_object(s, :struct_consts)
  s = load_object(s, :symbols)
  return lie_algebra(R, struct_consts, s; check=false)
end
