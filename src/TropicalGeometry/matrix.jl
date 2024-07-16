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
    is_tropically_generic(A::MatrixElem{<:TropicalSemiringElem}, minOrMax)

Check if a collection of vectors in the tropical torus (given as columns of a matrix `A`) are in tropical general position with respect to the `minOrMax` convention.

# Examples
```jldoctest
julia> A = matrix(tropical_semiring(),[1 0;0 1])
[(1)   (0)]
[(0)   (1)]

julia> is_tropically_generic(A,min)
true
```
"""
function is_tropically_generic(A::MatrixElem{<:TropicalSemiringElem}, minOrMax)
    @req convention(A) == minOrMax "Semiring convention not as declared"
    if ncols(A) == nrows(A)
        return Polymake.tropical.tsgn(A) != 0
    elseif ncols(A)>nrows(A)
	i = 0
        for b in subsets(ncols(A), nrows(A))
            if Polymake.tropical.tsgn(A[:,b]) == 0
                return false
		break
	    end
        end
	return i == 0
    else
	i = 0
        for b in subsets(nrows(A),ncols(A))
            if Polymake.tropical.tsgn(A[b,:]) == 0
                return false
		break
            end
        end
	return i == 0
    end
end
