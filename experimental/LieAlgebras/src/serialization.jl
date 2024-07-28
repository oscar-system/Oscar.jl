Oscar.encode_type(::Type{<: AbstractLieAlgebra}) = "AbstractLieAlgebra"
@register_serialization_type AbstractLieAlgebra uses_id

function save_object(s::SerializerState, L::AbstractLieAlgebra{T}) where T <: FieldElem
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_field)
    save_object(s, _struct_consts(L), :struct_const)
    save_object(s, symbols(L), :symbols)
  end
end
