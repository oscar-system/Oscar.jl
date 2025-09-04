#############################################################################
##
##  Bases of subspaces of a finite vector space

@doc raw"""
    bases_of_subspaces(V::AbstractAlgebra.Generic.FreeModule{T}, k::Int) where T <: FinFieldElem

Return an iterator of the matrices that are echelonized bases of the
`k`-dimensional subspaces of `V`.

# Examples
```jldoctest
julia> V = vector_space(GF(2), 3)
Vector space of dimension 3 over F

julia> for b in Oscar.bases_of_subspaces(V, 2) println(b); end
[1 0 0; 0 1 0]
[1 0 1; 0 1 0]
[1 0 0; 0 1 1]
[1 0 1; 0 1 1]
[1 0 0; 0 0 1]
[1 1 0; 0 0 1]
[0 1 0; 0 0 1]
```
"""
function bases_of_subspaces(V::AbstractAlgebra.Generic.FreeModule{T}, k::Int) where T <: FinFieldElem
  return SubspacesIterator{T}(V, k)
end

# For an iterator of `k`-dimensional subspaces in an `n`-dimensional space
# over the field with `q` elements,
# the state is `(comb_iter, choice, prod_iter, values)`
# where
# - `comb_iter` is the iterator `combinations(1:n, k)`,
# - `choice` is the current state of `comb_iter`,
#   a vector that denotes the `k` pivot column positions in the matrix,
# - `prod_iter` is an iterator of vectors, where the current vector
#   corresponds to those entries in the not necessarily zero positions
#   of the matrix that are not pivots, and
# - `values` is the current state ot `prod_iter`.

struct SubspacesIterator{T <: FinFieldElem}
  V::AbstractAlgebra.Generic.AbstractAlgebra.Generic.FreeModule{T}
  k::Int

  function SubspacesIterator{T}(V, k) where T
    @req 0 <= k <= vector_space_dim(V) "wrong dimension"
    return new(V, k)
  end
end

eltype(::SubspacesIterator{T}) where T = dense_matrix_type(T)

function length(S::SubspacesIterator{T}) where T
  n = vector_space_dim(S.V)
  if n < 2 * S.k
    k = n - S.k
  else
    k = S.k
  end
  q = order(base_ring(S.V))
  size = one(q)
  qn = q^n
  qd = q
  for i in 1:k
    size = div(size * (qn-1), (qd-1))
    qn = div(qn, q)
    qd = qd * q
  end
  return BigInt(size)
end

function _matrix_for_subspaces_iterator(F, k, n, choice, values)
  mat = zero_matrix(F, k, n)
  o = one(F)
  pos = 1
  for i in 1:k
    mat[i, choice[i]] = o
    for j in (choice[i]+1):n
      if !(j in choice)
        mat[i,j] = values[pos]
        pos = pos + 1
      end
    end
  end
  return mat
end

function Base.iterate(S::SubspacesIterator{T}) where T <: FinFieldElem
  F = base_ring(S.V)
  n = vector_space_dim(S.V)
  k = S.k

  comb = combinations(1:n, k)
  choice, stat_comb = iterate(comb)

  N = n*k - div(k*(k-1), 2) - sum(choice)
  proditer = Base.Iterators.ProductIterator(Tuple([F for i in 1:N]))
  prodval, prodstate = iterate(proditer)

  mat = _matrix_for_subspaces_iterator(F, k, n, choice, prodval)
  return mat, (comb, stat_comb, proditer, prodstate)
end

function Base.iterate(S::SubspacesIterator{T}, status::Tuple{Oscar.Combinations{UnitRange{Int}}, Vector{Int}, <:Base.Iterators.ProductIterator, <:Union{Tuple, Bool}}) where T <: FinFieldElem
  F = base_ring(S.V)
  n = vector_space_dim(S.V)
  k = S.k

  comb = status[1]
  stat_comb = status[2]
  choice = stat_comb
  proditer = status[3]
  prodstate = status[4]

  next = iterate(proditer, prodstate)
  if next === nothing
    next = iterate(comb, choice)
    next === nothing && return nothing
    choice, stat_comb = next

    N = n*k - div(k*(k-1), 2) - sum(choice)
    proditer = Base.Iterators.ProductIterator(Tuple([F for i in 1:N]))
    prodval, prodstate = iterate(proditer)
  else
    prodval, prodstate = next
  end

  mat = _matrix_for_subspaces_iterator(F, k, n, choice, prodval)
  return mat, (comb, stat_comb, proditer, prodstate)
end
