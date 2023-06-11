@doc raw"""
    morphism_from_cox_variety(variety::AbstractNormalToricVariety)

Return the quotient morphism from the Cox variety to the toric variety in question.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> morphism_from_cox_variety(F4)
A toric morphism
```
"""
@attr ToricMorphism function morphism_from_cox_variety(variety::AbstractNormalToricVariety)
    mapping_matrix = matrix(ZZ, rays(fan(variety)))
    max_cones_for_cox_variety = ray_indices(maximal_cones(variety))
    rays_for_cox_variety = matrix(ZZ, [[if i==j 1 else 0 end for j in 1:nrays(variety)] for i in 1:nrays(variety)])
    cox_variety = normal_toric_variety(polyhedral_fan(rays_for_cox_variety, max_cones_for_cox_variety))
    return toric_morphism(cox_variety, mapping_matrix, variety)
end


@doc raw"""
    cox_variety(variety::AbstractNormalToricVariety)

Return the Cox variety of the toric variety in question.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> cox_variety(F4)
Normal toric variety
```
"""
@attr AbstractNormalToricVariety function cox_variety(variety::AbstractNormalToricVariety)
    return domain(morphism_from_cox_variety(variety))
end
