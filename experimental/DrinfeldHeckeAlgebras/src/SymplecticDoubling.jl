################################################################################
# method for symplectic doubling of a matrix group over a field
#
# Cassandra Koenen, 2025
################################################################################

function symplectic_doubling(G::MatrixGroup{T}) where {T <: FieldElem}
  symplectic_generators = MatElem{T}[]

  for g in gens(G)
    A = matrix(g)
    h = block_diagonal_matrix([A, transpose(inv(A))])
    push!(symplectic_generators, h)
  end

  H = matrix_group(symplectic_generators)

  return H
end
