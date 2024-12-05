
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
    tensor_matrices_of_operators(L::LieAlgebra, highest_weight::WeightLatticeElem, operators::Vector{RootSpaceElem}) -> Vector{SMat{ZZRingElem}}

Calculates the action matrices of the operators in `operators` on
the tensor product of multiples of the fundamental modules (with multiplicities in `highest_weight`).
Note that the highest weight module with highest weight `highest_weight` is a submodule of this tensor product. 
We use multiples of fundamentals to reduce the total dimension of the ambient space
"""
function tensor_matrices_of_operators(
  L::LieAlgebra, highest_weight::WeightLatticeElem, operators::Vector{RootSpaceElem}
)
  R = root_system(L)
  mats = [zero_matrix(SMat, ZZ, 1) for _ in operators]
  for i in 1:rank(R)
    highest_weight_i = highest_weight[i]
    iszero(highest_weight_i) && continue
    mats = [
      _tensor_product(mat_temp, mat_wi) for (mat_temp, mat_wi) in zip(
        mats,
        matrices_of_operators(L, highest_weight_i * fundamental_weight(R, i), operators),
      )
    ]
  end
  return mats
end
