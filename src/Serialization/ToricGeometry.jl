################################################################################
# Toric varieties
@registerSerializationType(AffineNormalToricVariety)
@registerSerializationType(NormalToricVariety)

function save_internal(s::SerializerState, ntv::AbstractNormalToricVariety)
    return Dict(
        :polymakeNTV => save_type_dispatch(s, ntv.polymakeNTV)
    )
end

function load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where {T <: Union{NormalToricVariety, AffineNormalToricVariety}}
    pmobj = load_type_dispatch(s, Polymake.BigObject, dict[:polymakeNTV])
    return T(pmobj)
end


################################################################################
# Torus invariant divisors on toric varieties
@registerSerializationType(ToricDivisor)

function save_internal(s::SerializerState, td::ToricDivisor)
    return Dict(
        :toric_variety => save_type_dispatch(s, td.toric_variety),
        :coeffs => save_type_dispatch(s, td.coeffs),
    )
end

function load_internal(s::DeserializerState, ::Type{ToricDivisor}, dict::Dict)
    tv = load_unknown_type(s, dict[:toric_variety])
    coeffs = load_type_dispatch(s, Vector{fmpz}, dict[:coeffs])
    all = Polymake._lookup_multi(pm_object(tv), "DIVISOR")
    index = 0
    for i in 1:length(all)
        if Vector{fmpz}(all[i].COEFFICIENTS) == coeffs
            index = i
            break
        end
    end
    pmdiv = Polymake._lookup_multi(pm_object(tv), "DIVISOR", index-1)
    return ToricDivisor(pmdiv, tv, coeffs)
end
