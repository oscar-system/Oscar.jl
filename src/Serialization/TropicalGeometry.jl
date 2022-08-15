# Tropical Semiring

encodeType(::Type{<:TropicalSemiring}) = "TropicalSemiring"
reverseTypeMap["TropicalSemiring"] = TropicalSemiring

function save_internal(s::SerializerState, T::TropicalSemiring)
    println("hello")
    return Dict(
        :convention => save_type_dispatch(s, convention(T))
    )
end

function load_internal(s::DeserializerState, ::Type{TropicalSemiring}, dict::Dict)
    return nothing
end
