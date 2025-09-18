@testset "Algebraic Shifting" begin
  K = simplicial_complex([[1, 3] , [2, 3]])
  n = n_vertices(K)
  W = weyl_group(:A, n - 1)

  @testset "Symmetric Shifting" begin
    shift = Oscar.symmetric_shift(GF(2), K, longest_element(W))
    @test facets(shift) == [Set([2, 1]), Set([3, 2])]
  end
  
  @testset "Partial Shift Graph" begin
    s = gens(W)
    U = uniform_hypergraph(K, 2)
    all_shifts = partial_shift_graph_vertices(QQ, U, W)
    directed_graph, edge_labels = partial_shift_graph(QQ, all_shifts)
    @test collect(edges(directed_graph)) ==  [Edge(t...) for t in [[2, 1], [3, 1], [3, 2]]]
    @test issetequal(word.(edge_labels[2, 1]), word.([s[1], s[2] * s[1], s[1] * s[2] * s[1], s[1] * s[2]]))
  end
end
