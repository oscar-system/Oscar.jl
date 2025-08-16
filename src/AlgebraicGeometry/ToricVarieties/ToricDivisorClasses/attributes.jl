@doc raw"""
    divisor_class(tdc::ToricDivisorClass)

Return the element of the class group corresponding to the toric divisor class `tdc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> cl = class_group_with_map(P2)[1];

julia> tdc = toric_divisor_class(P2, cl([1]))
Divisor class on a normal toric variety

julia> divisor_class(tdc)
Abelian group element [1]
```
"""
divisor_class(tdc::ToricDivisorClass) = tdc.class


@doc raw"""
    toric_variety(tdc::ToricDivisorClass)

Return the toric variety on which the toric divisor class `tdc` is defined.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> cl = class_group_with_map(P2)[1];

julia> tdc = toric_divisor_class(P2, cl([1]))
Divisor class on a normal toric variety

julia> toric_variety(tdc)
Normal toric variety
```
"""
toric_variety(tdc::ToricDivisorClass) = tdc.toric_variety


@doc raw"""
    toric_divisor(tdc::ToricDivisorClass)

Construct a toric divisor corresponding to the toric divisor class  `tdc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> cl = class_group_with_map(P2)[1];

julia> tdc = toric_divisor_class(P2, cl([1]))
Divisor class on a normal toric variety

julia> toric_divisor(tdc)
Torus-invariant, prime divisor on a normal toric variety
```
"""
@attr ToricDivisor function toric_divisor(tdc::ToricDivisorClass)
    f = map_from_torusinvariant_weil_divisor_group_to_class_group(toric_variety(tdc))
    coeffs = vec([ZZRingElem(x) for x in preimage(f, divisor_class(tdc)).coeff])
    return toric_divisor(toric_variety(tdc), coeffs)
end
