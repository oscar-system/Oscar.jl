@doc Markdown.doc"""
    morphism_from_cox_variety(variety::AbstractNormalToricVariety)

This method returns the quotient morphism from the Cox variety to the toric variety in question.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> morphism_from_cox_variety(F4)
A toric morphism
```
"""
@attr ToricMorphism function morphism_from_cox_variety(variety::AbstractNormalToricVariety)
    mapping_matrix = matrix(ZZ, rays(fan(variety)))
    max_cones_for_cox_variety = ray_indices(maximal_cones(variety))
    rays_for_cox_variety = matrix(ZZ,[[if i==j 1 else 0 end for j in 1:nrays(variety)] for i in 1:nrays(variety)])
    cox_variety = NormalToricVariety(PolyhedralFan(rays_for_cox_variety, max_cones_for_cox_variety))
    return ToricMorphism(cox_variety, mapping_matrix, variety)
end
export morphism_from_cox_variety


@doc Markdown.doc"""
    cox_variety(variety::AbstractNormalToricVariety)

This method returns the Cox variety of the toric variety in question.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> cox_variety(F4)
A normal toric variety
```
"""
@attr AbstractNormalToricVariety function cox_variety(variety::AbstractNormalToricVariety)
    return domain(morphism_from_cox_variety(variety))
end
export cox_variety
