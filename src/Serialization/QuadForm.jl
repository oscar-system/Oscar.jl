############################################################
# QuadSpace
@register_serialization_type Hecke.QuadSpace uses_id

type_params(V::Hecke.QuadSpace) = TypeParams(Hecke.QuadSpace, parent(gram_matrix(V)))

function save_object(s::SerializerState, V::Hecke.QuadSpace)
  save_object(s, gram_matrix(V))
end

function load_object(s::DeserializerState, tp::TypeParams{<:Hecke.QuadSpace, <:MatSpace})
  params = parameters(tp)
  gram = load_object(s, TypeParams(MatElem, params))
  F = base_ring(params)
  return quadratic_space(F, gram; cached=false)
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

function load_object(s::DeserializerState, tp::TypeParams{ZZLat, <:Tuple{Vararg{Pair}}})
  mat_space = tp[:basis]
  B = load_object(s, TypeParams(elem_type(mat_space), mat_space))
  return lattice(tp[:ambient_space], B; check=false)
end

############################################################
# QuadSpaceWithIsom
@register_serialization_type QuadSpaceWithIsom                                                                                    
type_params(QS::QuadSpaceWithIsom) = TypeParams(
  QuadSpaceWithIsom,
  :quad_space => space(QS),
  :isom => parent(isometry(QS)),
  :order => TypeParams(typeof(order_of_isometry(QS)), nothing)
)

function save_object(s::SerializerState, QS::QuadSpaceWithIsom)
  save_data_dict(s) do
    save_object(s, isometry(QS), :isom)
    if !Base.issingletontype(typeof(order_of_isometry(QS)))
      save_object(s, order_of_isometry(QS), :order)
    end
  end
end

function load_object(s::DeserializerState, tp::TypeParams{QuadSpaceWithIsom, <:Tuple{Vararg{Pair}}})
  mat_space = tp[:isom]
  isom = load_object(s, TypeParams(elem_type(mat_space), mat_space), :isom)
  order_type = tp[:order]

  if Base.issingletontype(order_type)
    n = order_type()
  else
    n = load_object(s, Int, :order)
  end
  return QuadSpaceWithIsom(tp[:quad_space], isom, n)
end

############################################################
# ZZLatWithIsom
@register_serialization_type ZZLatWithIsom

type_params(x::ZZLatWithIsom) = TypeParams(
  ZZLatWithIsom,
  :ambient_space => ambient_space(x),
  :basis => parent(basis_matrix(x))
)

function save_object(s::SerializerState, L::ZZLatWithIsom)
  save_object(s, basis_matrix(L))
end

function load_object(s::DeserializerState, tp::TypeParams{ZZLatWithIsom, <:Tuple{Vararg{Pair}}})
  mat_space = tp[:basis]
  B = load_object(s, TypeParams(elem_type(mat_space), mat_space))
  return lattice(tp[:ambient_space], B; check=false)
end
