# This file implements some basics to work with the Hermitian inner product over
# number fields.

# Ulrich Thiel, 2024



###########################################################################################
# Complex conjugation
###########################################################################################
# In Hecke v0.30.2 complex_conjugation was extended to more fields following my suggestion.
# I will only do one minor extension following Tommy Hofmanns suggestion:
complex_conjugation(K::QQAbField) = conj
complex_conjugation(K::QQBarField) = conj


###########################################################################################
# Hermitian scalar product
###########################################################################################
function scalar_product(v::AbstractAlgebra.Generic.FreeModuleElem{T}, w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

  @req parent(v) === parent(w) "Incompatible vector spaces"

  V = parent(v)
  K = base_ring(V)
  n = dim(V)
  s = zero(K)
  conj = complex_conjugation(K)
  for i=1:n
    s += v[i]*conj(w[i])
  end
  return s
end

###########################################################################################
# Check if a matrix is unitary
###########################################################################################
function is_orthogonal(M::MatElem)

  return is_square(M) && is_one(M*transpose(M))

end

function is_unitary(M::QQMatrix)

  return is_orthogonal(M)

end

function is_unitary(M::MatElem{T}) where T <: QQAlgFieldElem

  if !is_square(M)
    return false
  end

  # create the conjugate transpose of M
  K = base_ring(M)
  n = ncols(M)

  conj = complex_conjugation(K)
  Mct = transpose(map_entries(conj, M))

  return is_one(M*Mct)

end

function _is_nondegenerate_hermitian_matrix(J::MatElem{T}) where T <: QQAlgFieldElem
  is_square(J) || return false
  conjugation = complex_conjugation(base_ring(J))
  conjugate_transpose = transpose(map_entries(conjugation, J))
  return J == conjugate_transpose && !is_zero(det(J))
end

"""
    is_unitary(M::MatElem)
    is_unitary(G::MatGroup)
    is_unitary(M::MatElem, J::MatElem)
    is_unitary(G::MatGroup, J::MatElem)

Return whether `M`, or every generator of `G`, preserves the standard
Hermitian form or, when `J` is supplied, the Hermitian form with Gram matrix
`J`.

Matrices and vectors use OSCAR's right-action convention, so invariance means
`M * J * conjugate_transpose(M) == J`, where conjugation is the exact complex
conjugation of the coefficient field. The matrix `J` must be Hermitian and
nondegenerate. This predicate does not certify positivity at an archimedean
embedding; when `J` is known to be positive definite, it tests unitarity in the
usual inner-product sense.
"""
function is_unitary(M::MatElem{T}, J::MatElem{T}) where T <: QQAlgFieldElem
  if !is_square(M) ||
      nrows(M) != nrows(J) ||
      base_ring(M) !== base_ring(J) ||
      !_is_nondegenerate_hermitian_matrix(J)
    return false
  end

  conjugation = complex_conjugation(base_ring(M))
  conjugate_transpose = transpose(map_entries(conjugation, M))
  return M * J * conjugate_transpose == J
end

function is_unitary(M::MatGroupElem{T}) where T <: QQAlgFieldElem

  return is_unitary(matrix(M))

end

function is_unitary(G::MatGroup{T}) where T <: QQAlgFieldElem

  return all(is_unitary, gens(G))

end


function is_unitary(G::MatGroup{T}, J::MatElem{T}) where T <: QQAlgFieldElem
  if degree(G) != nrows(J) ||
      base_ring(G) !== base_ring(J) ||
      !_is_nondegenerate_hermitian_matrix(J)
    return false
  end
  return all(g -> is_unitary(matrix(g), J), gens(G))
end
