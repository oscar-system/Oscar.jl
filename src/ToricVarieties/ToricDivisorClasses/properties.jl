@attr Bool is_trivial(tdc::ToricDivisorClass) = iszero(divisor_class(tdc))
export is_trivial


@doc Markdown.doc"""
    is_effective(tdc::ToricDivisorClass)

Determines whether the toric divisor class `tdc` is effective, that is if a toric divisor
in this divisor class is linearly equivalent to an effective toric divisor.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety,2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, [1])
A divisor class on a normal toric variety

julia> is_effective(tdc)
true

julia> tdc2 = ToricDivisorClass(P2, [-1])
A divisor class on a normal toric variety

julia> is_effective(tdc2)
false
```
"""
@attr Bool is_effective(tdc::ToricDivisorClass) = is_effective(toric_divisor(tdc))
export is_effective
