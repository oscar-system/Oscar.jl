"""
    Oscar.LieAlgebras.coefficient_vector(M::MatElem{T}, basis::Vector{<:MatElem{T}}) where {T}

Return the vector of coefficients of the matrix `M` in the basis `basis`.
Requires that `basis` is linearly independent and that `M` lies in the span of `basis`.

# Examples
```jldoctest
julia> Oscar.LieAlgebras.coefficient_vector(matrix(QQ, [1 2;3 4]), [matrix(QQ, [1 0;0 0]), matrix(QQ, [0 1;0 2]), matrix(QQ, [0 0;-1 0])])
[1   2   -3]
```
"""
function coefficient_vector(M::MatElem{T}, basis::Vector{<:MatElem{T}}) where {T<:FieldElem}
  (nr, nc) = size(M)
  all(b -> size(b) == (nr, nc), basis) ||
    throw(DimensionMismatch("The matrices in B must have the same size as M."))
  @req all(b -> base_ring(b) == base_ring(M), basis) "Incompatible base rings"
  lgs = zero(M, length(basis), nr * nc)
  for (k, b) in enumerate(basis)
    for i in 1:nr
      lgs[k, (i - 1) * nc .+ (1:nc)] = view(b, i:i, :)
    end
  end
  rhs = similar(M, 1, nr * nc)
  for i in 1:nr
    rhs[1, (i - 1) * nc .+ (1:nc)] = view(M, i:i, :)
  end
  sol = solve(lgs, rhs; side=:left)
  return sol
end
