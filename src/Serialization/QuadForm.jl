############################################################
# QuadSpace
@register_serialization_type Hecke.QuadSpace uses_params

type_params(V::Hecke.QuadSpace) = type_params(gram_matrix(V))

function save_object(s::SerializerState, V::Hecke.QuadSpace)
  save_object(s, gram_matrix(V))
end

function load_object(s::DeserializerState, ::Type{<:Hecke.QuadSpace}, params::MatSpace)
  gram = load_object(s, MatElem, params)
  F =  base_ring(params)
  return quadratic_space(F, gram)
end

############################################################
# ZZLat
@register_serialization_type ZZLat

function save_object(s::SerializerState, L::ZZLat)
  save_data_dict(s) do
    save_typed_object(s, basis_matrix(L), :basis)
    save_typed_object(s, ambient_space(L), :ambient_space)
  end
end

function load_object(s::DeserializerState, ::Type{ZZLat})
  B = load_typed_object(s, :basis)
  V = load_typed_object(s, :ambient_space)
  return lattice(V, B)
end
