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
function det(A::Generic.MatSpaceElem{<:TropicalSemiringElem})
    detA = zero(base_ring(A))
    nrows(A)!=ncols(A) && return detA # return tropical zero if matrix not square
    for sigma in AbstractAlgebra.SymmetricGroup(nrows(A)) # otherwise follow Leibniz Formula
        detA += prod([ A[i,sigma[i]] for i in 1:nrows(A) ])
    end
    return detA
end

function det(A::Matrix{TropicalSemiringElem{minOrMax}}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    return det(matrix(TropicalSemiring{minOrMax},A))
end

@doc raw"""
    t_generic(A::Generic.MatSpaceElem{<: TropicalSemiringElem}, maxormin)
Checks if a collection of vectors in the tropical torus (given as columns of a matrix `A`) are in tropical general position with respect to the `maxormin` convention.
# Examples
```jldoctest
julia> A = matrix(tropcial_semiring(),[1 0;01])
[(1)   (0)]
[(0)   (1)]

julia> t_generic(A,min)
true
```
"""
function tropically_generic(A::Generic.MatSpaceElem{<:TropicalSemiringElem}, minOrMax)
    @req convention(A) == minOrMax "Semiring convention not as declared"
    if ncols(A) == nrows(A)
        return Polymake.tropical.tsgn(A) != 0
    elseif ncols(A)>nrows(A)
        for b in subsets(ncols(A), nrows(A))
            if Polymake.tropical.tsgn(A[:,b]) == 0
                return false
            end
        end
    else
        for b in subsets(nrows(A),ncols(A))
            if Polymake.tropical.tsgn(A[b,:]) == 0
                return false
            end
        end
    end
end
