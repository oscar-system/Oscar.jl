# Tropical Semiring

encodeType(::Type{TropicalSemiring{S}}) where S = "TropicalSemiring{$S}"
reverseTypeMap["TropicalSemiring{typeof(min)}"] = TropicalSemiring{typeof(min)}
reverseTypeMap["TropicalSemiring{typeof(max)}"] = TropicalSemiring{typeof(max)}

## elements
encodeType(::Type{<:TropicalSemiringElem})= "TropicalSemiringElem"
reverseTypeMap["TropicalSemiringElem"] = TropicalSemiringElem

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
encodeType(::Type{<: TropicalHypersurface}) = "TropicalHypersurface"
reverseTypeMap["TropicalHypersurface"] = TropicalHypersurface

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
encodeType(::Type{<: TropicalCurve}) = "TropicalCurve"
reverseTypeMap["TropicalCurve"] = TropicalCurve

function save_internal(s::SerializerState, t_curve::TropicalCurve{M, EMB}) where {M, EMB}
  if EMB
    return Dict(
        :polyhedral_complex => save_type_dispatch(s, underlying_polyhedral_complex(t_curve))
    )
  else
    return Dict(
        :graph => save_type_dispatch(s, graph(t_curve))
    )
  end
end

function load_internal(s::DeserializerState,
                       ::Type{<: TropicalCurve},
                       dict::Dict)
  if haskey(dict, :polyhedral_complex)
    return TropicalCurve(
      load_type_dispatch(s, PolyhedralComplex{fmpq}, dict[:polyhedral_complex])
    )
  else
    return TropicalCurve(
      load_type_dispatch(s, IncidenceMatrix, dict[:graph])
    )
  end
end
