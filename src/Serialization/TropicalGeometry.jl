# Tropical Semiring
function save_internal(s::SerializerState, T::TropicalSemiring)
    return Dict(
        :convention => save_type_dispatch()
    )
end

function load_internal(s::DeserializerState, ::Type{TropicalSemiring}, dict::Dict)
    return nothing
end
