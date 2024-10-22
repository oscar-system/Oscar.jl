@testset "Algebraic Shifting" begin
  K = simplicial_complex([[1, 2, 3], [2, 4]])

  @testset "Partial Shift Graph" begin
    n = n_vertices(K)
    W = weyl_group(:A, n - 1)
    partial_shift_graph_nodes(QQ, K, collect(W))
  end
end
