@doc raw"""
    morphism_from_cox_variety(variety::NormalToricVarietyType)

Return the quotient morphism from the Cox variety to the toric variety in question.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> morphism_from_cox_variety(F4)
Toric morphism
```
"""
@attr ToricMorphism function morphism_from_cox_variety(variety::NormalToricVarietyType)
    mapping_matrix = matrix(ZZ, rays(variety))
    max_cones_for_cox_variety = ray_indices(maximal_cones(variety))
    rays_for_cox_variety = matrix(ZZ, [[if i==j 1 else 0 end for j in 1:n_rays(variety)] for i in 1:n_rays(variety)])
    cox_variety = normal_toric_variety(polyhedral_fan(max_cones_for_cox_variety, rays_for_cox_variety))
    return toric_morphism(cox_variety, mapping_matrix, variety)
end


@doc raw"""
    cox_variety(variety::NormalToricVarietyType)

Return the Cox variety of the toric variety in question.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4)
Normal toric variety

julia> cox_variety(F4)
Normal toric variety
```
"""
@attr NormalToricVarietyType function cox_variety(variety::NormalToricVarietyType)
    return domain(morphism_from_cox_variety(variety))
end
