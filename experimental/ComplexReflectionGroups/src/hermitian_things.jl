# This file implements some basics to work with the Hermitian inner product over
# number fields.

# Ulrich Thiel, 2024



###########################################################################################
# Complex conjugation
###########################################################################################
# In Hecke v0.30.2 complex_conjugation was extended to more fields following my suggestion.
# I will only do one minor extension following Tommy Hofmanns suggestion:
complex_conjugation(K::QQAbField) = conj


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

function is_unitary(M::MatElem{T}) where T <: NumFieldElem

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

function is_unitary(M::MatrixGroupElem{T}) where T <: QQAlgFieldElem

  return is_unitary(matrix(M))

end

function is_unitary(G::MatrixGroup{T}) where T <: QQAlgFieldElem
  
  return all(is_unitary, gens(G))

end
