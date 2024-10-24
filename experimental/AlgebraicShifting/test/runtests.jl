@testset "Algebraic Shifting" begin
  K = simplicial_complex([[1, 3] , [2, 3]])

  @testset "Partial Shift Graph" begin
    n = n_vertices(K)
    W = weyl_group(:A, n - 1)
    s = gens(W)
    all_shifts = partial_shift_graph_vertices(QQ, K, W)
    directed_graph, edge_labels = partial_shift_graph(QQ, vertices)
    @test collect(edges(directed_graph)) ==  [Edge(t...) for t in [[2, 1], [3, 1], [3, 2]]]
    @test word.(edge_labels[2, 1]) == word.([s[1], s[1] * s[2], s[1] * s[2] * s[1], s[2] * s[1]])
  end
end
