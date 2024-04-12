@testset "GaussianGraphicalModels" begin

  S = gaussian_ring(3)
  G = graph_from_edges(Directed, [[1,2],[2,3]])
  M = graphical_model(G, S)
  
end



@testset "DiscreteGraphicalModels" begin
  
end



@testset "ConditionalIndependence" begin
  
end



@testset "Phylogenetics" begin
  
end





