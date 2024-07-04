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

Return the tropical determinant of `A`.

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
    detA = zero(base_ring(A))
    nrows(A)!=ncols(A) && return detA # return tropical zero if matrix not square
    for sigma in AbstractAlgebra.SymmetricGroup(nrows(A)) # otherwise follow Leibniz Formula
        detA += prod([ A[i,sigma[i]] for i in 1:nrows(A) ])
    end
    return detA
end

function det(A::Matrix{R}) where {R<:Union{TropicalSemiringElem,MPolyRingElem{<:TropicalSemiringElem},PolyRingElem{<:TropicalSemiringElem}}}
    return det(matrix(parent(first(A)),A))
end
