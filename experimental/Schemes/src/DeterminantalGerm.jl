################################################################################
## Some more getter functions
################################################################################

function defining_matrix(X::DeterminantalGerm)
  return X.A
end

function determinantal_type(X::DeterminantalGerm)
  n, m = size(X.A)
  return n, m, X.t
end

################################################################################
## More constructors
################################################################################

#TODO: To be added

################################################################################
## T1_GL module
################################################################################

function _R_ij(A::MatElem, i::Integer, j::Integer)
  R_ij = zero(A)
  R_ij[i, :] = A[j, :]
  return R_ij
end

function _C_ij(A::MatElem, i::Integer, j::Integer)
  C_ij = zero(A)
  C_ij[:, i] = A[:, j]
  return C_ij
end

@attr SubquoModule function T1_GL_module(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}})
  # Use the transpose for _vec() below to match the entries correctly
  A = transpose(defining_matrix(X))
  m, n = size(A)
  L = base_ring(parent(A))
  N = ngens(L)
  modulus = vcat([derivative.(A, i) for i in 1:N],
                 [_C_ij(A, i, j) for i in 1:n for j in 1:n],
                 [_R_ij(A, i, j) for i in 1:m for j in 1:m]
                )
  # index switched, since (m, n) is the size of the transpose
  F = FreeMod(L, [Symbol("E[$i,$j]") for i in 1:n for j in 1:m])
  return quo(F, F.(_vec.(modulus)))[1]
end

tjurina_GL_number(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}}) = vector_space_dim(T1_GL_module(X))

is_determinantally_rigid(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}}) = is_zero(T1_GL_module(X))

@attr Bool is_EIDS(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}}) = krull_dim(T1_GL_module(X)) <= 0

basis_versal_det_unfolding(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}}) = vector_space_basis(T1_GL_module(X))

