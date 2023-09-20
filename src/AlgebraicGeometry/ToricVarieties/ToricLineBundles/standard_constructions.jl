########################
# 1: Special attributes of toric varieties
########################

@doc raw"""
    structure_sheaf(v::NormalToricVarietyType)

Construct the structure sheaf of a normal toric variety.
For convenience, we also support `structure_sheaf(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> structure_sheaf(v)
Toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle structure_sheaf(v::NormalToricVarietyType) = toric_line_bundle(v, zero(picard_group(v)))


@doc raw"""
    anticanonical_bundle(v::NormalToricVarietyType)

Construct the anticanonical bundle of a normal toric variety.
For convenience, we also support `anticanonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> anticanonical_bundle(v)
Toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle anticanonical_bundle(v::NormalToricVarietyType) = toric_line_bundle(v, anticanonical_divisor(v))


@doc raw"""
    canonical_bundle(v::NormalToricVarietyType)

Construct the canonical bundle of a normal toric variety.
For convenience, we also support `canonical_bundle(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> canonical_bundle(v)
Toric line bundle on a normal toric variety
```
"""
@attr ToricLineBundle canonical_bundle(v::NormalToricVarietyType) = inv(anticanonical_bundle(v))
