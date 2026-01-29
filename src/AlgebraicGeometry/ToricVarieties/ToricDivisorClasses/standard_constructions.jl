########################
# 1: Special attributes of toric varieties
########################

@doc raw"""
    trivial_divisor_class(v::NormalToricVarietyType)

Construct the trivial divisor class of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> trivial_divisor_class(v)
Divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass trivial_divisor_class(v::NormalToricVarietyType) =
  toric_divisor_class(
    trivial_divisor(v)
  )

@doc raw"""
    anticanonical_divisor_class(v::NormalToricVarietyType)

Construct the anticanonical divisor class of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> anticanonical_divisor_class(v)
Divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass anticanonical_divisor_class(v::NormalToricVarietyType) =
  -canonical_divisor_class(v)

@doc raw"""
    canonical_divisor_class(v::NormalToricVarietyType)

Construct the canonical divisor class of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> canonical_divisor_class(v)
Divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass canonical_divisor_class(v::NormalToricVarietyType) =
  toric_divisor_class(
    canonical_divisor(v)
  )
