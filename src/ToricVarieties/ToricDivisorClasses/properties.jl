@doc Markdown.doc"""
    istrivial(tdc::ToricDivisorClass)

Checks if the divisor class `td` is trivial.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, [1])
A divisor class on a normal toric variety

julia> istrivial(tdc)
false
```
"""
function istrivial(tdc::ToricDivisorClass)
    return get_attribute!(tdc, :istrivial) do
        return iszero(divisor_class(tdc))::Bool
    end
end
export istrivial
