@testset "Graphical Models tests" begin
  tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3]])
  model = jukes_cantor_model(tree)
  @test is_isomorphic(graph(model), graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
  @test Oscar.number_states(model) == 4
  @test model isa Oscar.PhylogeneticModel
end
