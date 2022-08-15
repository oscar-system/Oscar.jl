# Tropical Semiring

encodeType(::Type{TropicalSemiring{S}}) where S = "TropicalSemiring{$S}"
reverseTypeMap["TropicalSemiring{typeof(min)}"] = TropicalSemiring{typeof(min)}
reverseTypeMap["TropicalSemiring{typeof(max)}"] = TropicalSemiring{typeof(max)}

# elements
function save_internal(s::SerializerState, t::TropicalSemiringElem{S}) where S
    T = parent(t)
    return Dict(
        :parent => save_type_dispatch(s, T),
        
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.NmodRing}, dict::Dict)
    modulus = load_type_dispatch(s, UInt64, dict[:modulus])
    return Nemo.NmodRing(modulus)
end


