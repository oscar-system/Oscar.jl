########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    StructureSheaf(v::AbstractNormalToricVariety)

Construct the structure sheaf of a normal toric variety.
For convenience, we also support `structure_sheaf(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> StructureSheaf(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle function StructureSheaf(v::AbstractNormalToricVariety)
    class = zero(picard_group(v))
    return ToricLineBundle(v, class)
end
structure_sheaf(v::AbstractNormalToricVariety) = StructureSheaf(v)
export StructureSheaf, structure_sheaf


@doc Markdown.doc"""
    AnticanonicalBundle(v::AbstractNormalToricVariety)

Construct the anticanonical bundle of a normal toric variety.
For convenience, we also support `anticanonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> AnticanonicalBundle(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle function AnticanonicalBundle(v::AbstractNormalToricVariety)
    return ToricLineBundle(v, sum(cox_ring(v).d))
end
anticanonical_bundle(v::AbstractNormalToricVariety) = AnticanonicalBundle(v)
export AnticanonicalBundle, anticanonical_bundle


@doc Markdown.doc"""
    CanonicalBundle(v::AbstractNormalToricVariety)

Construct the canonical bundle of a normal toric variety.
For convenience, we also support `canonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> CanonicalBundle(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle function CanonicalBundle(v::AbstractNormalToricVariety)
    return ToricLineBundle(v, (-1)*sum(cox_ring(v).d))
end
canonical_bundle(v::AbstractNormalToricVariety) = CanonicalBundle(v)
export CanonicalBundle, canonical_bundle
