########################
# 1: Special attributes of toric varieties
########################

@doc raw"""
    trivial_divisor(v::AbstractNormalToricVariety)

Construct the trivial divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> trivial_divisor(v)
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor trivial_divisor(v::AbstractNormalToricVariety) = toric_divisor(v, zeros(ZZRingElem, nrays(v)))


@doc raw"""
    anticanonical_divisor(v::AbstractNormalToricVariety)

Construct the anticanonical divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> anticanonical_divisor(v)
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor anticanonical_divisor(v::AbstractNormalToricVariety) = toric_divisor(v, fill(ZZRingElem(1), nrays(v)))


@doc raw"""
    canonical_divisor(v::AbstractNormalToricVariety)

Construct the canonical divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> canonical_divisor(v)
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor canonical_divisor(v::AbstractNormalToricVariety) = toric_divisor(v, fill(ZZRingElem(-1), nrays(v)))
