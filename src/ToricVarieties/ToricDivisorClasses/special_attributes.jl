########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    trivial_divisor_class(v::AbstractNormalToricVariety)

Construct the trivial divisor class of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> trivial_divisor_class(v)
A divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass trivial_divisor_class(v::AbstractNormalToricVariety) = ToricDivisorClass(trivial_divisor(v))
export trivial_divisor_class


@doc Markdown.doc"""
    anticanonical_divisor_class(v::AbstractNormalToricVariety)

Construct the anticanonical divisor class of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> anticanonical_divisor_class(v)
A divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass anticanonical_divisor_class(v::AbstractNormalToricVariety) = ToricDivisorClass(anticanonical_divisor(v))
export anticanonical_divisor_class


@doc Markdown.doc"""
    canonical_divisor_class(v::AbstractNormalToricVariety)

Construct the canonical divisor class of a normal toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> canonical_divisor_class(v)
A divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass canonical_divisor_class(v::AbstractNormalToricVariety) = ToricDivisorClass(canonical_divisor(v))
export canonical_divisor_class
