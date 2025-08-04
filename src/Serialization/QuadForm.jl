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
  return lattice(params[:ambient_space], B; check=false)
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

function load_object(s::DeserializerState, ::Type{QuadSpaceWithIsom}, params::Dict)
  quad_space = params[:quad_space]
  mat_space = params[:isom]
  isom = load_object(s, MatElem{elem_type(mat_space)}, mat_space, :isom)
  order_type = params[:order]

  if Base.issingletontype(order_type)
    n = order_type()
  else
    n = load_object(s, Int, :order)
  end
  return QuadSpaceWithIsom(quad_space, isom, n)
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

function load_object(s::DeserializerState, ::Type{ZZLatWithIsom}, params::Dict)
  quad_space = params[:ambient_space]
  mat_space = params[:basis]
  B = load_object(s, MatElem{elem_type(mat_space)}, mat_space)
  return lattice(quad_space, B; check=false)
end
