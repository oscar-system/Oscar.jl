@attr Bool is_trivial(tdc::ToricDivisorClass) = iszero(divisor_class(tdc))


@doc raw"""
    is_effective(tdc::ToricDivisorClass)

Determines whether the toric divisor class `tdc` is effective, that is if a toric divisor
in this divisor class is linearly equivalent to an effective toric divisor.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety,2)
Normal toric variety

julia> tdc = toric_divisor_class(P2, [1])
Divisor class on a normal toric variety

julia> is_effective(tdc)
true

julia> tdc2 = toric_divisor_class(P2, [-1])
Divisor class on a normal toric variety

julia> is_effective(tdc2)
false
```
"""
@attr Bool is_effective(tdc::ToricDivisorClass) = (length(basis_of_global_sections(toric_line_bundle(tdc))) > 0)
