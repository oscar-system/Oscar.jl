########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    trivial_divisor(v::AbstractNormalToricVariety)

Construct the trivial divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> trivial_divisor(v)
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor trivial_divisor(v::AbstractNormalToricVariety) = ToricDivisor(v, zeros(fmpz, nrays(v)))
export trivial_divisor


@doc Markdown.doc"""
    anticanonical_divisor(v::AbstractNormalToricVariety)

Construct the anticanonical divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> anticanonical_divisor(v)
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor anticanonical_divisor(v::AbstractNormalToricVariety) = ToricDivisor(v, fill(fmpz(1), nrays(v)))
export anticanonical_divisor


@doc Markdown.doc"""
    canonical_divisor(v::AbstractNormalToricVariety)

Construct the canonical divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> canonical_divisor(v)
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor canonical_divisor(v::AbstractNormalToricVariety) = ToricDivisor(v, fill(fmpz(-1), nrays(v)))
export canonical_divisor
