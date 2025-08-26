# This file implements explicit models of symplectic reflection groups.
#
# At the moment there's nothing really implemented yet.
#
# Ulrich Thiel, 2024


function symplectic_doubling(G::MatrixGroup)

  K = base_ring(G)
  n = degree(G)
  matspace = matrix_space(K, 2*n, 2*n)
  gens_symp = elem_type(matspace)[]

  for g in gens(G)
    g_symp = matspace(block_diagonal_matrix([matrix(g), transpose(matrix(g^-1))]))
    push!(gens_symp, g_symp)
  end

  G_symp = matrix_group(gens_symp)

  set_attribute!(G_symp, :order, order(G))

  return G_symp

end

function symplectic_reflection_group(W::MatrixGroup)

  if !is_complex_reflection_group(W)
    throw(ArgumentError("Group must be a complex reflection group"))
  end

  W_symp = symplectic_doubling(W)

  set_attribute!(W_symp, :is_symplectic_reflection_group, true)

  return W_symp
end

function is_symplectic_reflection_group(G::MatrixGroup)

  if has_attribute(G, :is_symplectic_reflection_group)
    return get_attribute(G, :is_symplectic_reflection_group)
  end
  
  #this should be upgraded later to work with a general matrix group
  error("Not implemented yet")

end
