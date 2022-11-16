########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    structure_sheaf(v::AbstractNormalToricVariety)

Construct the structure sheaf of a normal toric variety.
For convenience, we also support `structure_sheaf(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> structure_sheaf(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle structure_sheaf(v::AbstractNormalToricVariety) = ToricLineBundle(v, zero(picard_group(v)))
export structure_sheaf


@doc Markdown.doc"""
    anticanonical_bundle(v::AbstractNormalToricVariety)

Construct the anticanonical bundle of a normal toric variety.
For convenience, we also support `anticanonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> anticanonical_bundle(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle anticanonical_bundle(v::AbstractNormalToricVariety) = prod(ToricLineBundle(v, d) for d in torusinvariant_prime_divisors(v))
export anticanonical_bundle


@doc Markdown.doc"""
    canonical_bundle(v::AbstractNormalToricVariety)

Construct the canonical bundle of a normal toric variety.
For convenience, we also support `canonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> canonical_bundle(v)
A toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle canonical_bundle(v::AbstractNormalToricVariety) = inv(anticanonical_bundle(v))
export canonical_bundle
