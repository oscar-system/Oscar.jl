# Tropical Semiring
@registerSerializationType(TropicalSemiring{typeof(min)})
@registerSerializationType(TropicalSemiring{typeof(max)})
has_elem_basic_encoding(obj::TropicalSemiring) = true

## elements
@registerSerializationType(TropicalSemiringElem)
type_needs_params(T::Type{<: TropicalSemiringElem}) = true

function save_type_params(s::SerializerState, x::TropicalSemiringElem, key::Symbol)
  s.key = key
  data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_typed_object(s, parent(x), :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:TropicalSemiringElem}, dict::Dict)
  return load_typed_object(s, dict)
end

function load_type_params(s::DeserializerState, ::Type{<:TropicalSemiringElem},
                          parent::TropicalSemiring)
  return load_typed_object(s, dict)
end

function load_object_with_params(s::DeserializerState,
                                 ::Type{<:TropicalSemiringElem},
                                 str::String, params::TropicalSemiring)
  return params(load_object(s, QQFieldElem, str))
end

# Tropical Hypersurfaces
@registerSerializationType(TropicalHypersurface, true)

function save_internal(s::SerializerState, t_surf::TropicalHypersurface)
  return Dict(
    :tropical_polynomial => save_type_dispatch(s, polynomial(t_surf))
  )
end

function load_internal(s::DeserializerState,
                       ::Type{<: TropicalHypersurface},
                       dict::Dict)
  polynomial = load_type_dispatch(s, MPolyRingElem, dict[:tropical_polynomial])
  return TropicalHypersurface(polynomial)
end

# Tropical Curves
@registerSerializationType(TropicalCurve, true)

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
      load_type_dispatch(s, PolyhedralComplex{QQFieldElem}, dict[:polyhedral_complex])
    )
  else
    return TropicalCurve(
      load_type_dispatch(s, IncidenceMatrix, dict[:graph])
    )
  end
end
