@testset "GaussianGraphicalModels" begin

  S = gaussian_ring(3)
  s = gens(S)
  G = graph_from_edges(Directed, [[1,2],[2,3]])
  M = graphical_model(G, S)

  @test vanishing_ideal(M) == ideal([-s[1, 2]*s[2, 3] + s[1, 3]*s[2, 2]])
end






