########################
# 1: Special attributes of toric varieties
########################

@doc raw"""
    trivial_divisor(v::NormalToricVarietyType)

Construct the trivial divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> trivial_divisor(v)
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor trivial_divisor(v::NormalToricVarietyType) = toric_divisor(
  v, zeros(ZZRingElem, n_rays(v))
)

@doc raw"""
    anticanonical_divisor(v::NormalToricVarietyType)

Construct the anticanonical divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> anticanonical_divisor(v)
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor anticanonical_divisor(v::NormalToricVarietyType) = -canonical_divisor(v)

@doc raw"""
    canonical_divisor(v::NormalToricVarietyType)

Construct the canonical divisor of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> canonical_divisor(v)
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
@attr ToricDivisor canonical_divisor(v::NormalToricVarietyType) = toric_divisor(
  v, fill(ZZRingElem(-1), n_rays(v))
)
