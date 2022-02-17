@doc Markdown.doc"""
    divisor_class(tdc::ToricDivisorClass)

Return the element of the class group corresponding to the toric divisor class `tdc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, class_group(P2)([1]))
A divisor class on a normal toric variety

julia> divisor_class(tdc)
Element of
GrpAb: Z
with components [1]
```
"""
divisor_class(tdc::ToricDivisorClass) = tdc.class
export divisor_class


@doc Markdown.doc"""
    toric_variety(tdc::ToricDivisorClass)

Return the toric variety on which the toric divisor class `tdc` is defined.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, class_group(P2)([1]))
A divisor class on a normal toric variety

julia> toric_variety(tdc)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
toric_variety(tdc::ToricDivisorClass) = tdc.toric_variety
export toric_variety


@doc Markdown.doc"""
    toric_divisor(tdc::ToricDivisorClass)

Constructs a toric divisor corresponding to the toric divisor class  `tdc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, class_group(P2)([1]))
A divisor class on a normal toric variety

julia> toric_divisor(tdc)
A torus-invariant, prime divisor on a normal toric variety
```
"""
@attr ToricDivisor function toric_divisor(tdc::ToricDivisorClass)
    f = map_from_torusinvariant_weil_divisor_group_to_class_group(toric_variety(tdc))
    coeffs = vec([fmpz(x) for x in preimage(f, divisor_class(tdc)).coeff])
    return ToricDivisor(toric_variety(tdc), coeffs)
end
export toric_divisor
