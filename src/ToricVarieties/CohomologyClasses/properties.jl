@doc Markdown.doc"""
    istrivial(c::CohomologyClass)

Checks if the cohomology class is trivial.

# Examples
```jldoctest
julia> F5 = hirzebruch_surface(5)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> cc = CohomologyClass(ToricDivisorClass(F5,[1,2]))
A cohomology class on a normal toric variety given by -3*t2 + x2

julia> istrivial(cc)
false
```
"""
@attr Bool istrivial(cc::CohomologyClass) = iszero(polynomial(cc))
export istrivial
