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
function divisor_class(tdc::ToricDivisorClass)
    return tdc.class
end
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
function toric_variety(tdc::ToricDivisorClass)
    return tdc.toric_variety
end
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
function toric_divisor(tdc::ToricDivisorClass)
    return get_attribute!(tdc, :toric_divisor) do
        f = map_from_torus_invariant_weil_divisor_group_to_class_group(toric_variety(tdc))
        coeffs = vec([fmpz(x) for x in preimage(f, divisor_class(tdc)).coeff])
        return ToricDivisor(toric_variety(tdc), coeffs)
    end
end
export toric_divisor
