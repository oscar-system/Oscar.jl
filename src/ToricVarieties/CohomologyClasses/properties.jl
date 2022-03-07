@doc Markdown.doc"""
    istrivial(c::CohomologyClass)

Checks if the cohomology class is trivial.

# Examples
```jldoctest
julia> H5 = hirzebruch_surface(5)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> cc = CohomologyClass(ToricDivisorClass(H5,[1,2]))
A cohomology class on a normal toric variety given by -9*x3 + 2*x4

julia> istrivial(cc)
false
```
"""
@attr Bool istrivial(cc::CohomologyClass) = iszero(polynomial(cc))
export istrivial
