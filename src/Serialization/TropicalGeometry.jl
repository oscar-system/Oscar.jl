# Tropical Semiring

encodeType(::Type{TropicalSemiring{S}}) where S = "TropicalSemiring{$S}"
reverseTypeMap["TropicalSemiring{typeof(min)}"] = TropicalSemiring{typeof(min)}
reverseTypeMap["TropicalSemiring{typeof(max)}"] = TropicalSemiring{typeof(max)}

# elements
encodeType(::Type{TropicalSemiringElem{S}}) where S = "TropicalSemiringElem{$S}"
reverseTypeMap["TropicalSemiringElem{typeof(min)}"] = TropicalSemiringElem{typeof(min)}
reverseTypeMap["TropicalSemiringElem{typeof(max)}"] = TropicalSemiringElem{typeof(max)}

function save_internal(s::SerializerState, t::TropicalSemiringElem{S}) where {S}
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

function load_internal_internal(s::DeserializerState,
                       ::Type{TropicalSemiringElem{S}},
                       dict::Dict) where S
  parent = load_type_dispatch(s, TropicalSemiring{S}, dict[:parent])
  return parent(load_type_dispatch(s, fmpq, dict[:data]))
end




