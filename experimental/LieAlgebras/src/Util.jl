"""
    coefficient_vector(M::MatElem{T}, basis::Vector{<:MatElem{T}}) where {T}

Return the vector of coefficients of the matrix `M` in the basis `basis`.
Requires that `basis` is linearly independent and that `M` lies in the span of `basis`.

# Examples
```jldoctest; setup = :(using Oscar.LieAlgebras)
julia> coefficient_vector(matrix(QQ, [1 2;3 4]), [matrix(QQ, [1 0;0 0]), matrix(QQ, [0 1;0 2]), matrix(QQ, [0 0;-1 0])])
[1   2   -3]
```
"""
function coefficient_vector(M::MatElem{T}, basis::Vector{<:MatElem{T}}) where {T<:FieldElem}
  (nr, nc) = size(M)
  all(b -> size(b) == (nr, nc), basis) ||
    throw(DimensionMismatch("The matrices in B must have the same size as M."))
  lgs = similar(M, nr * nc, length(basis))
  for (k, b) in enumerate(basis)
    for i in 1:nr, j in 1:nc
      lgs[(i - 1) * nc + j, k] = b[i, j]
    end
  end
  rhs = similar(M, nr * nc, 1)
  for i in 1:nr, j in 1:nc
    rhs[(i - 1) * nc + j, 1] = M[i, j]
  end
  sol = solve(lgs, rhs; side = :right)
  return transpose(sol)
end
