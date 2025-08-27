function matrices_of_operators(
  L::LieAlgebra, highest_weight::WeightLatticeElem, operators::Vector{RootSpaceElem}
)
  # used in tensor_matrices_of_operators
  R = root_system(L)
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L, highest_weight)
  dimV = size(struct_consts, 2)
  transformation_matrices = [
    sparse_matrix(coefficient_ring(L), dimV, dimV) for _ in 1:number_of_roots(R)
  ]
  for i in 1:number_of_roots(R), j in 1:dimV
    transformation_matrices[i][j] = struct_consts[i, j]
  end

  matrices_of_operators = map(operators) do op
    is_pos_root, pos_i = is_positive_root_with_index(op)
    is_neg_root, neg_i = is_negative_root_with_index(op)
    if is_pos_root
      return change_base_ring(ZZ, transformation_matrices[pos_i + number_of_positive_roots(R)])  # take f_alpha for positive root alpha
    else if is_neg_root
      return change_base_ring(ZZ, transformation_matrices[neg_i - number_of_positive_roots(R)])  # take e_alpha for negative root alpha
    else
      error("unreachable")
    end
  end

  return matrices_of_operators
end
