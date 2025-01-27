# Tropical Semiring
@register_serialization_type TropicalSemiring{typeof(min)}
@register_serialization_type TropicalSemiring{typeof(max)}

## elements
@register_serialization_type TropicalSemiringElem uses_params

function save_object(s::SerializerState, x::TropicalSemiringElem)
  str = string(x)
  save_data_basic(s, String(strip(str, ['(', ')'])))
end

function load_object(s::DeserializerState, ::Type{<:TropicalSemiringElem}, params::Vector)
  t_ring = params[end]
  load_node(s) do str
    if str == "∞" || str == "-∞" || str == "infty" || str == "-infty"
      return inf(t_ring)
    else
      return t_ring(load_object(s, QQFieldElem))
    end
  end
end

function load_object(s::DeserializerState, T::Type{<:TropicalSemiringElem},
                     parent_ring::TropicalSemiring)
  return load_object(s, T, get_parents(parent_ring))
end

# Tropical Hypersurfaces
@register_serialization_type TropicalHypersurface uses_id uses_params

type_params(t::T) where T <: TropicalHypersurface = type_params(tropical_polynomial(t))

function save_object(s::SerializerState, t::T) where T <: TropicalHypersurface
  save_data_dict(s) do
    save_object(s, tropical_polynomial(t), :tropical_polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{<: TropicalHypersurface},
                     params::MPolyRing)
  polynomial = load_object(s, MPolyRingElem, params, :tropical_polynomial)
  return tropical_hypersurface(polynomial)
end

# Tropical Curves
@register_serialization_type TropicalCurve uses_id uses_params

type_params(t::TropicalCurve{M, true}) where M = type_params(polyhedral_complex(t))
# here to handle annnoying edge case
type_params(t::TropicalCurve{M, false}) where M = "graph"

function save_object(s::SerializerState, t::TropicalCurve{M, EMB}) where {M, EMB}
  save_data_dict(s) do
    if EMB
      save_object(s, polyhedral_complex(t), :polyhedral_complex)
    else
      save_object(s, graph(t), :graph)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: TropicalCurve}, params::Field)
  return tropical_curve(
    load_object(s, PolyhedralComplex, params, :polyhedral_complex)
  )
end

function load_object(s::DeserializerState, ::Type{<: TropicalCurve}, ::String)
  return tropical_curve(
    load_object(s, Graph{Undirected}, :graph)
  )
end
