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

function load_object(s::DeserializerState,
                                 ::Type{<:TropicalSemiringElem},
                                 str::String, params::TropicalSemiring)
  return params(load_object(s, QQFieldElem, str))
end

# Tropical Hypersurfaces
@registerSerializationType(TropicalHypersurface, true)

function save_object(s::SerializerState, t_surf::T) where T <: TropicalHypersurface
  save_typed_object(s, polynomial(t_surf), :data)
end

function load_object(s::DeserializerState, ::Type{<: TropicalHypersurface}, dict::Dict)
  polynomial = load_typed_object(s, dict)
  return TropicalHypersurface(polynomial)
end

# Tropical Curves
@registerSerializationType(TropicalCurve, true)

function save_object(s::SerializerState, t_curve::TropicalCurve{M, EMB}) where {M, EMB}
  data_dict(s) do 
    if EMB
      save_typed_object(s, underlying_polyhedral_complex(t_curve), :polyhedral_complex)
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
    return TropicalCurve(
      load_typed_object(s, dict[:polyhedral_complex])
    )
  else
    return TropicalCurve(
      load_typed_object(s, dict[:graph])
    )
  end
end
