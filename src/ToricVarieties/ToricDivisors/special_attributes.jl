########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    TrivialDivisor(v::AbstractNormalToricVariety)

Construct the trivial divisor of a normal toric variety.
For convenience, we also support `trivial_divisor(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> TrivialDivisor(v)
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor function TrivialDivisor(v::AbstractNormalToricVariety)
    return ToricDivisor(v, zeros(fmpz, nrays(v)))
end
trivial_divisor(v::AbstractNormalToricVariety) = TrivialDivisor(v)
export TrivialDivisor, trivial_divisor


@doc Markdown.doc"""
    AnticanonicalDivisor(v::AbstractNormalToricVariety)

Construct the anticanonical divisor of a normal toric variety.
For convenience, we also support `anticanonical_divisor(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> AnticanonicalDivisor(v)
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor function AnticanonicalDivisor(v::AbstractNormalToricVariety)
    return ToricDivisor(v, [fmpz(1) for i in 1:nrays(v)])
end
anticanonical_divisor(v::AbstractNormalToricVariety) = AnticanonicalDivisor(v)
export AnticanonicalDivisor, anticanonical_divisor


@doc Markdown.doc"""
    CanonicalDivisor(v::AbstractNormalToricVariety)

Construct the canonical divisor of a normal toric variety.
For convenience, we also support `canonical_divisor(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> CanonicalDivisor(v)
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor function CanonicalDivisor(v::AbstractNormalToricVariety)
    return ToricDivisor(v, [fmpz(-1) for i in 1:nrays(v)])
end
canonical_divisor(v::AbstractNormalToricVariety) = CanonicalDivisor(v)
export CanonicalDivisor, canonical_divisor
