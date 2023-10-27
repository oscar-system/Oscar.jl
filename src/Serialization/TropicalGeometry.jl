# Tropical Semiring
@register_serialization_type TropicalSemiring{typeof(min)}
@register_serialization_type TropicalSemiring{typeof(max)}

## elements
@register_serialization_type TropicalSemiringElem uses_params

function save_type_params(s::SerializerState, x::T) where {T <: TropicalSemiringElem}
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_typed_object(s, parent(x), :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:TropicalSemiringElem}, dict::Dict{Symbol, Any})
  return load_typed_object(s, dict)
end

function load_type_params(s::DeserializerState, ::Type{<:TropicalSemiringElem},
                          parent::TropicalSemiring)
  return load_typed_object(s, dict)
end

function load_object(s::DeserializerState,
                                 ::Type{<:TropicalSemiringElem},
                                 str::String, params::TropicalSemiring)
  if str == "∞" || str == "-∞" || str == "infty" || str == "-infty"
    return inf(params)
  else
    # looks like (q)
    return params(load_object(s, QQFieldElem, String(strip(str, ['(', ')']))))
  end
end

# Tropical Hypersurfaces
@register_serialization_type TropicalHypersurface uses_id

function save_object(s::SerializerState, t_surf::T) where T <: TropicalHypersurface
  save_data_dict(s) do
    save_typed_object(s, tropical_polynomial(t_surf), :tropical_polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{<: TropicalHypersurface}, dict::Dict)
  polynomial = load_typed_object(s, dict[:tropical_polynomial])
  return tropical_hypersurface(polynomial)
end

# Tropical Curves
@register_serialization_type TropicalCurve uses_id

function save_object(s::SerializerState, t_curve::TropicalCurve{M, EMB}) where {M, EMB}
  save_data_dict(s) do
    if EMB
      save_typed_object(s, polyhedral_complex(t_curve), :polyhedral_complex)
      save_object(s, true, :is_embedded)
    else
      save_typed_object(s, graph(t_curve), :graph)
      save_object(s, false, :is_embedded)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: TropicalCurve}, dict::Dict)
  EMB = parse(Bool, dict[:is_embedded])
  if EMB
    return tropical_curve(
      load_typed_object(s, dict[:polyhedral_complex])
    )
  else
    return tropical_curve(
      load_typed_object(s, dict[:graph])
    )
  end
end
