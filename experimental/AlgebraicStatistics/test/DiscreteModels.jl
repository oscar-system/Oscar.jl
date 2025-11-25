@testset "DiscreteGraphicalModels" begin
  mktempdir() do path
    @testset "Directed" begin
      DG = graph_from_edges(Directed, [[1,2], [1,3]])
      M = discrete_graphical_model(DG, [2,2,2])
      V1 = vanishing_ideal(M; algorithm=:f4)
      @test iszero(V1)

      test_save_load_roundtrip(path, M) do loaded
        @test loaded == M
      end
    end

    @testset "Undirected" begin
      UG = graph_from_edges([[1, 2], [1, 3]])
      M = discrete_graphical_model(UG, [2,2,2])
      V2 = vanishing_ideal(M)
      _, p = model_ring(M)
      @test gens(V2) == [p[2,1,1] * p[2,2,2] - p[2,2,1] * p[2,1,2],  p[1,1,1]*p[1,2,2] - p[1,2,1]*p[1,1,2]]

      test_save_load_roundtrip(path, M) do loaded
        @test loaded == M
      end
    end
  end
end
