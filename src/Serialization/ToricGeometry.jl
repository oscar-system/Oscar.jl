function save_internal(s::SerializerState, ntv::AbstractNormalToricVariety)
    return Dict(
        :polymakeNTV => save_type_dispatch(s, ntv.polymakeNTV)
    )
end

function load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where {T <: Union{NormalToricVariety, AffineNormalToricVariety}}
    pmobj = load_type_dispatch(s, Polymake.BigObject, dict[:polymakeNTV])
    return T(pmobj)
end
