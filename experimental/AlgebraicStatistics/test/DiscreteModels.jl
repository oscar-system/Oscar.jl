@testset "DiscreteGraphicalModels" begin
  M = discrete_graphical_model(graph_from_edges([[1,2], [2,3]]), [2,2,2])
  parameter_ring(M)
  model_ring(M)

  # still need directed discrete case
end
