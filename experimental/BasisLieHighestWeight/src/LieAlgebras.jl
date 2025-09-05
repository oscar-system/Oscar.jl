function matrices_of_operators(
  L::LieAlgebra, highest_weight::WeightLatticeElem, operators::Vector{RootSpaceElem}
)
  # used in tensor_matrices_of_operators
  R = root_system(L)
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L, highest_weight)
  dimV = size(struct_consts, 2)
  transformation_matrices_e = [
    sparse_matrix(coefficient_ring(L), dimV, dimV) for _ in 1:number_of_positive_roots(R)
  ]
  transformation_matrices_f = [
    sparse_matrix(coefficient_ring(L), dimV, dimV) for _ in 1:number_of_positive_roots(R)
  ]
  for i in 1:number_of_positive_roots(R), j in 1:dimV
    transformation_matrices_e[i][j] = struct_consts[i, j]
    transformation_matrices_f[i][j] = struct_consts[i + number_of_positive_roots(R), j]
  end

  matrices_of_operators = map(operators) do op
    # operators act downwards, so we need to swap e and f
    if ((fl, i) = is_positive_root_with_index(op); fl)
      return change_base_ring(ZZ, transformation_matrices_f[i])
    elseif ((fl, i) = is_negative_root_with_index(op); fl)
      return change_base_ring(ZZ, transformation_matrices_e[i])
    else
      error("unreachable")
    end
  end

  return matrices_of_operators
end
