################################################
# 1: Intersections involving closed subvarieties
################################################

@doc Markdown.doc"""
    RationalEquivalenceClass(v::AbstractNormalToricVariety, coefficients::Vector{T}) where {T <: IntegerUnion}

Construct the rational equivalence class of algebraic cycles corresponding to a linear combination of cones.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(P2, [1,2,3])
A torus-invariant, non-prime divisor on a normal toric variety

julia> RationalEquivalenceClass(d)
A rational equivalence class on a normal toric variety represented by 6V(x3)
```
"""

function Base.:*(ac::RationalEquivalenceClass, sv::ClosedSubvarietyOfToricVariety)
    if toric_variety(ac) !== toric_variety(sv)
        throw(ArgumentError("The rational equivalence class and the closed subvariety must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ac * RationalEquivalenceClass(sv)
end


function Base.:*(sv::ClosedSubvarietyOfToricVariety, ac::RationalEquivalenceClass)
    if toric_variety(ac) !== toric_variety(sv)
        throw(ArgumentError("The rational equivalence class and the closed subvariety must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ac * RationalEquivalenceClass(sv)
end
