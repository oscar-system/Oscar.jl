############################################################
# QuadSpace
@register_serialization_type Hecke.QuadSpace uses_id

type_params(V::Hecke.QuadSpace) = TypeParams(Hecke.QuadSpace, parent(gram_matrix(V)))

function save_object(s::SerializerState, V::Hecke.QuadSpace)
  save_object(s, gram_matrix(V))
end

function load_object(s::DeserializerState, ::Type{<:Hecke.QuadSpace}, params::MatSpace)
  gram = load_object(s, MatElem, params)
  F = base_ring(params)
  return quadratic_space(F, gram)
end

############################################################
# ZZLat
@register_serialization_type ZZLat

type_params(L::ZZLat) = TypeParams(
  ZZLat,
  :basis => parent(basis_matrix(L)),
  :ambient_space => ambient_space(L)
)

function save_object(s::SerializerState, L::ZZLat)
  save_object(s, basis_matrix(L))
end

function load_object(s::DeserializerState, ::Type{ZZLat}, params::Dict)
  mat_space = params[:basis]
  B = load_object(s, elem_type(mat_space), mat_space)
  return lattice(params[:ambient_space], B)
end
