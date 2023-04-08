"""
    coefficient_vector(M::MatElem{T}, basis::Vector{<:MatElem{T}}) where {T}

Returns the vector of coefficients of the matrix `M` in the basis `basis`.
Requires that `basis` is linearly independent and that `M` lies in the span of `basis`.

# Examples
```jldoctest; setup = :(using Oscar.LieAlgebras)
julia> coefficient_vector(matrix(QQ, [1 2;3 4]), [matrix(QQ, [1 0;0 0]), matrix(QQ, [0 1;0 2]), matrix(QQ, [0 0;-1 0])])
[1   2   -3]
```
"""
function coefficient_vector(
  M::MatElem{T}, basis::Vector{<:MatElem{T}}
) where {T<:RingElement}
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
  return transpose(solve(lgs, rhs))
end

"""
    is_valid_dynkin(dynkin::Char, n::Int)

Returns true, if there given parameters uniquely define a dynkin diagram,
i.e. are of one of the forms
  * ``A_n`` for ``n \\geq 1``,
  * ``B_n`` for ``n \\geq 2``,
  * ``C_n`` for ``n \\geq 2``,
  * ``D_n`` for ``n \\geq 4``,
  * ``E_5``, ``E_6``, ``E_7``,
  * ``F_4``,
  * ``G_2``.
"""
function is_valid_dynkin(dynkin::Char, n::Int)
  if dynkin == 'A'
    return n >= 1
  elseif dynkin == 'B'
    return n >= 2
  elseif dynkin == 'C'
    return n >= 2
  elseif dynkin == 'D'
    return n >= 4
  elseif dynkin == 'E'
    return 6 <= n <= 8
  elseif dynkin == 'F'
    return n == 4
  elseif dynkin == 'G'
    return n == 2
  else
    return false
  end
end
