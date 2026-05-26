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

# TODO: move to a more fitting place (maybe AbstractAlgebra)
function _sym_mat_gens(A::MatElem)
  n, m = size(A)
  @req n == m "Matrix 'A' must be a quadratic."
  R = base_ring(A)
  gens = typeof(A)[]
  sizehint!(gens, n*(n+1)/2)
  for i in 1:n
    for j in i:n
      tmp = zero(A)
      tmp[i,j] = one(R)
      tmp[j,i] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end

# TODO: move to a more fitting place (maybe AbstractAlgebra)
function _skew_sym_mat_gens(A::MatElem)
  n, m = size(A)
  @req n == m "Matrix 'A' must be a quadratic."
  R = base_ring(A)
  characteristic(R) == 2 && return _sym_mat_gens(A)
  gens = typeof(A)[]
  sizehint!(gens, n*(n-1)/2)
  for i in 1:n
    for j in i+1:n
      tmp = zero(A)
      tmp[i,j] = -one(R)
      tmp[j,i] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end


@attr SubquoModule function T1_GL_module(X::DeterminantalGerm)
  # transposing, since '_vec' vcats the columms of A and we would rather read rowwise
  A = transpose(defining_matrix(X))
  m, n = size(A)
  L = base_ring(parent(A))
  N = ngens(L)
  # defining_matrix has size n x m
  F = FreeMod(L, [Symbol("E[$i,$j]") for i in 1:n for j in 1:m])

  if X.type === Val{:generic}()
    rels = vcat([derivative.(A, i) for i in 1:N],
                [_C_ij(A, i, j) for i in 1:n for j in 1:n],
                [_R_ij(A, i, j) for i in 1:m for j in 1:m]
               )
    return SubquoModule(F, gens(F), F.(_vec.(rels)))
  end

  # Case: X.type === Val{:symmetric} or # X.type === Val{:skew_symmetric}
  erz = X.type === Val{:symmetric}() ? _sym_mat_gens(A) : _skew_sym_mat_gens(A)
  rels = vcat([derivative.(A, i) for i in 1:N],
              [_R_ij(A, i, j) + _C_ij(A, i, j) for i in 1:n for j in 1:n],
             )
  return SubquoModule(F, F.(_vec.(erz)), F.(_vec.(rels)))
end

tjurina_GL_number(X::DeterminantalGerm) = vector_space_dim(T1_GL_module(X))

is_determinantally_rigid(X::DeterminantalGerm) = is_zero(T1_GL_module(X))

@attr Bool is_EIDS(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}}) = krull_dim(T1_GL_module(X)) <= 0

basis_versal_det_unfolding(X::DeterminantalGerm) = vector_space_basis(T1_GL_module(X))

