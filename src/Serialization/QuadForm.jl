############################################################
# QuadSpace
@register_serialization_type Hecke.QuadSpace uses_id

type_and_params(V::Hecke.QuadSpace) = TypeAndParams(Hecke.QuadSpace, parent(gram_matrix(V)))

function save_object(s::SerializerState, V::Hecke.QuadSpace)
  save_object(s, gram_matrix(V))
end

function load_object(s::DeserializerState, tp::TypeAndParams{<:Hecke.QuadSpace, <:MatSpace})
  params = params(tp)
  gram = load_object(s, TypeAndParams(MatElem, params))
  F = base_ring(params)
  return quadratic_space(F, gram; cached=false)
end

############################################################
# ZZLat
@register_serialization_type ZZLat

type_and_params(L::ZZLat) = TypeAndParams(
  ZZLat,
  :basis => parent(basis_matrix(L)),
  :ambient_space => ambient_space(L)
)

function save_object(s::SerializerState, L::ZZLat)
  save_object(s, basis_matrix(L))
end

function load_object(s::DeserializerState, tp::TypeAndParams{ZZLat, <:Tuple{Vararg{Pair}}})
  mat_space = tp[:basis]
  B = load_object(s, TypeAndParams(elem_type(mat_space), mat_space))
  return lattice(tp[:ambient_space], B; check=false)
end

############################################################
# QuadSpaceWithIsom
@register_serialization_type QuadSpaceWithIsom                                                                                    
type_and_params(QS::QuadSpaceWithIsom) = TypeAndParams(
  QuadSpaceWithIsom,
  :quad_space => space(QS),
  :isom => parent(isometry(QS)),
  :order => type_and_params(order_of_isometry(QS))
)

function save_object(s::SerializerState, QS::QuadSpaceWithIsom)
  save_data_dict(s) do
    save_object(s, isometry(QS), :isom)
    if !Base.issingletontype(typeof(order_of_isometry(QS)))
      save_object(s, order_of_isometry(QS), :order)
    end
  end
end

function load_object(s::DeserializerState, tp::TypeAndParams{QuadSpaceWithIsom, <:Tuple{Vararg{Pair}}})
  mat_space = tp[:isom]
  isom = load_object(s, TypeAndParams(elem_type(mat_space), mat_space), :isom)
  order_type = type(tp[:order])
  if Base.issingletontype(order_type)
    n = order_type()
  else
    n = load_object(s, tp[:order], :order)
  end
  return QuadSpaceWithIsom(tp[:quad_space], isom, n)
end

############################################################
# ZZLatWithIsom
@register_serialization_type ZZLatWithIsom

type_and_params(x::ZZLatWithIsom) = TypeAndParams(
  ZZLatWithIsom,
  :ambient_space => ambient_space(x),
  :basis => parent(basis_matrix(x))
)

function save_object(s::SerializerState, L::ZZLatWithIsom)
  save_object(s, basis_matrix(L))
end

function load_object(s::DeserializerState, tp::TypeAndParams{ZZLatWithIsom, <:Tuple{Vararg{Pair}}})
  mat_space = tp[:basis]
  B = load_object(s, TypeAndParams(elem_type(mat_space), mat_space))
  return lattice(tp[:ambient_space], B; check=false)
end
