########################
# 1: Special attributes of toric varieties
########################

@doc Markdown.doc"""
    TrivialDivisorClass(v::AbstractNormalToricVariety)

Construct the trivial divisor class of a normal toric variety.
For convenience, we also support `trivial_divisor_class(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> TrivialDivisorClass(v)
A divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass function TrivialDivisorClass(v::AbstractNormalToricVariety)
    class = zero(class_group(v))
    return ToricDivisorClass(v, class)
end
trivial_divisor_class(v::AbstractNormalToricVariety) = TrivialDivisorClass(v)
export TrivialDivisorClass, trivial_divisor_class


@doc Markdown.doc"""
    AnticanonicalDivisorClass(v::AbstractNormalToricVariety)

Construct the anticanonical divisor class of a normal toric variety.
For convenience, we also support `anticanonical_divisor_class(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> AnticanonicalDivisorClass(v)
A divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass function AnticanonicalDivisorClass(v::AbstractNormalToricVariety)
    return ToricDivisorClass(v, sum(cox_ring(v).d))
end
anticanonical_divisor_class(v::AbstractNormalToricVariety) = AnticanonicalDivisorClass(v)
export AnticanonicalDivisorClass, anticanonical_divisor_class


@doc Markdown.doc"""
    CanonicalDivisorClass(v::AbstractNormalToricVariety)

Construct the canonical divisor class of a normal toric variety.
For convenience, we also support `canonical_divisor_class(variety)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> CanonicalDivisorClass(v)
A divisor class on a normal toric variety
```
"""
@attr ToricDivisorClass function CanonicalDivisorClass(v::AbstractNormalToricVariety)
    return ToricDivisorClass(v, (-1)*sum(cox_ring(v).d))
end
canonical_divisor_class(v::AbstractNormalToricVariety) = CanonicalDivisorClass(v)
export CanonicalDivisorClass, canonical_divisor_class
