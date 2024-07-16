################################################################################
#
#  Tropical matrices
#
################################################################################


################################################################################
#
#  Substraction-less tropical alternatives to generic functions
#
################################################################################

@doc raw"""
    det(A::Generic.MatSpaceElem{<: TropicalSemiringElem})

Return the tropical determinant of `A`.  That is, this function evaluates the tropicalization of the ordinary determinant considered as a multivariate polynomial.

That computation is equivalent to solving a linear assignment problem from combinatorial optimization.  The implementation employs the Hungarian method, which is polynomial time.  See Chapter 3 in [Jos21](@cite).

!!! note
    This function effectively overwrites the `det` command for tropical matrices.  This means that functions like `minors` will use the tropical determinant when used on a tropical matrix.

# Examples
```jldoctest
julia> A = matrix(tropical_semiring(),[1 2; 3 4])
[(1)   (2)]
[(3)   (4)]

julia> det(A)
(5)
```
"""
function det(A::Generic.MatSpaceElem{R}) where {R<:Union{TropicalSemiringElem,MPolyRingElem{<:TropicalSemiringElem},PolyRingElem{<:TropicalSemiringElem}}}
    T = base_ring(A)
    nrows(A)!=ncols(A) && return zero(T) # return tropical zero if matrix not square
    return T(Polymake.tropical.tdet(A))
end

function det(A::Matrix{R}) where {R<:Union{TropicalSemiringElem,MPolyRingElem{<:TropicalSemiringElem},PolyRingElem{<:TropicalSemiringElem}}}
    return det(matrix(parent(first(A)),A))
end
