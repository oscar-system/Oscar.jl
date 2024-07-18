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
    det(A::MatrixElem{<: TropicalSemiringElem})

Return the tropical determinant of `A`.  That is, this function evaluates the tropicalization of the ordinary determinant considered as a multivariate polynomial at `A`.

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
function det(A::MatrixElem{<:TropicalSemiringElem})
  @req nrows(A) == ncols(A) "Non-square matrix"
  T = base_ring(A)
  return T(Polymake.tropical.tdet(A))
end

function det(A::Matrix{<:TropicalSemiringElem})
  @req 0 < nrows(A) == ncols(A) "Non-square or empty matrix"
  return det(matrix(parent(first(A)),A))
end
