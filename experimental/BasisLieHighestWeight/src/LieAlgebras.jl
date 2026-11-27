function matrices_of_operators(
  L::LieAlgebra, highest_weight::WeightLatticeElem, operators::Vector{RootSpaceElem}
)
  # used in tensor_matrices_of_operators
  R = root_system(L)
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L, highest_weight)
  dimV = size(struct_consts, 2)

  return map(operators) do op
    # operators act downwards, so we need to swap e and f
    if ((fl, i) = is_positive_root_with_index(op); fl)
      basis_index = number_of_positive_roots(R) + i # f_i
    elseif ((fl, i) = is_negative_root_with_index(op); fl)
      basis_index = i # e_i
    else
      error("unreachable")
    end
    mat = sparse_matrix(ZZ, dimV, dimV)
    for j in 1:dimV
      mat[j] = change_base_ring(ZZ, struct_consts[basis_index, j])
    end
    return mat
  end
end
