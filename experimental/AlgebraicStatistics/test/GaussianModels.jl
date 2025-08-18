const ColoredGGM{Directed} = GaussianGraphicalModel{
  Directed,
  @NamedTuple{color::S}
} where S <: Oscar.GraphMap;
# 

@testset "GaussianGraphicalModels" begin
  DG = graph_from_edges(Directed, [[1,2],[2,3]])
  @testset "Directed" begin
    M1 = gaussian_graphical_model(DG)
    cov_mat = covariance_matrix(M1)
    V1 = vanishing_ideal(M1)
    @test V1 == ideal(
      [-cov_mat[1, 2] * cov_mat[2, 3] + cov_mat[1, 3] * cov_mat[2, 2]]
    )

    label!(DG,
           Dict((1, 2) => "pink", (2, 3) => "pink"),
           Dict(i => "green" for i in 1:3);
           name=:color)
    M2 = gaussian_graphical_model(DG)
    cov_mat = covariance_matrix(M2)
    V2 = vanishing_ideal(M2)
    @test V2 == ideal(
      [-cov_mat[1, 2] * cov_mat[2, 3] + cov_mat[1, 3] * cov_mat[2, 2]]
    )

  end

  UG = complete_bipartite_graph(1, 2)
  @testset "Undirected" begin
    M3 = gaussian_graphical_model(UG)
    #TODO fix this when algorithm = :eliminate ?
    V3 = vanishing_ideal(M3; algorithm=:kernel)
    cov_mat = covariance_matrix(M3)
    @test V3 == ideal([
      -cov_mat[1, 1] * cov_mat[2, 3] + cov_mat[1, 2] * cov_mat[1, 3]])
  end

  @testset "Mixed" begin
    MG = mixed_graph_from_edges(collect(edges(DG)), collect(edges(UG)))
    M4 = gaussian_graphical_model(MG)
  end
end
