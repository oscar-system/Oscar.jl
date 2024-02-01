
function _tensor_product(A, B)
  return kronecker_product(A, identity_matrix(SMat, ZZ, nrows(B))) +
         kronecker_product(identity_matrix(SMat, ZZ, nrows(A)), B)
end

function _tensor_power(A, k)
  return sum(
    kronecker_product(
      kronecker_product(identity_matrix(SMat, ZZ, nrows(A)^(j - 1)), A),
      identity_matrix(SMat, ZZ, nrows(A)^(k - j)),
    ) for j in 1:k
  )
end

@doc raw"""
    tensor_matrices_of_operators(L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}, operators::Vector{GAP.Obj}) -> Vector{SMat{ZZRingElem}}

Calculates the action matrices of the operators in `operators` on
the tensor product of multiples of the fundamental modules (with multiplicities in `highest_weight`).
Note that the highest weight module with highest weight `highest_weight` is a submodule of this tensor product. 
We use multiples of fundamentals to reduce the total dimension of the ambient space
"""
function tensor_matrices_of_operators(
  L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}, operators::Vector{GAP.Obj}
)
  matrices_of_operators = [zero_matrix(SMat, ZZ, 1) for _ in operators]
  for (i, highest_weight_i) in enumerate(Int.(highest_weight))
    if highest_weight_i <= 0
      continue
    end
    wi = ZZ.(1:length(highest_weight) .== i) # i-th fundamental weight
    matrices_of_operators = [
      _tensor_product(mat_temp, mat_wi) for (mat_temp, mat_wi) in zip(
        matrices_of_operators,
        matrices_of_operators_gap(L, highest_weight_i * wi, operators),
      )
    ]
  end
  return matrices_of_operators
end
