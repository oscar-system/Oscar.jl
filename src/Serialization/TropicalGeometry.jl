# Tropical Semiring
@registerSerializationType(TropicalSemiring{typeof(min)})
@registerSerializationType(TropicalSemiring{typeof(max)})

## elements
@registerSerializationType(TropicalSemiringElem)

function save_internal(s::SerializerState, t::TropicalSemiringElem)
  T = parent(t)
  return Dict(
    :parent => save_type_dispatch(s, T),
    :data => save_type_dispatch(s, data(t))
  )
end

function load_internal(s::DeserializerState,
                       ::Type{TropicalSemiringElem{S}},
                       dict::Dict) where S
  parent = load_type_dispatch(s, TropicalSemiring{S}, dict[:parent])
  return parent(load_type_dispatch(s, fmpq, dict[:data]))
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{TropicalSemiringElem{S}},
                                   dict::Dict,
                                   parent::TropicalSemiring{S}) where S
  return parent(load_type_dispatch(s, fmpq, dict[:data]))
end


# Tropical Hypersurfaces
@registerSerializationType(TropicalHypersurface)

function save_internal(s::SerializerState, t_surf::TropicalHypersurface)
    return Dict(
        :tropical_polynomial => save_type_dispatch(s, polynomial(t_surf))
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: TropicalHypersurface},
                       dict::Dict)
  polynomial = load_type_dispatch(s, MPolyElem, dict[:tropical_polynomial])
  return TropicalHypersurface(polynomial)
end

# Tropical Curves
@registerSerializationType(TropicalCurve)

function save_internal(s::SerializerState, t_curve::TropicalCurve{M, EMB}) where {M, EMB}
  if EMB
    return Dict(
        :polyhedral_complex => save_type_dispatch(s, underlying_polyhedral_complex(t_curve)),
        :is_embedded => save_type_dispatch(s, true)
    )
  else
    return Dict(
        :graph => save_type_dispatch(s, graph(t_curve)),
        :is_embedded => save_type_dispatch(s, false)
    )
  end
end

function load_internal(s::DeserializerState,
                       ::Type{<: TropicalCurve},
                       dict::Dict)
  EMB = load_type_dispatch(s, Bool, dict[:is_embedded])
  if EMB
    return TropicalCurve(
      load_type_dispatch(s, PolyhedralComplex{fmpq}, dict[:polyhedral_complex])
    )
  else
    return TropicalCurve(
      load_type_dispatch(s, IncidenceMatrix, dict[:graph])
    )
  end
end
