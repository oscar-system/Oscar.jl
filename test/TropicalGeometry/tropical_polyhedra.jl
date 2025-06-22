
@testset "TropicalPolyhedron{min}" begin
  TT = tropical_semiring(min)
  pts1_a = TT[0 -4 -1; 0 -1 0]
  P1_a = tropical_convex_hull(pts1_a)
  @test P1_a isa TropicalPolyhedron{typeof(min)}
  @test dim(P1_a) == 1

  ## Adding redundant points gives same tropical polytope with same properties
  pts1_b = TT[0 -4 -1; 0 -1 0; 0 -3 0]
  P1_b = tropical_convex_hull(pts1_b)
  @test issetequal(P1_a |> vertices, P1_b |> vertices)
  @test dim(P1_b) == 1

  cov1 = maximal_covectors(P1_a)
  @test issetequal(cov1, [IncidenceMatrix([1], [2], [2]), IncidenceMatrix([1], [2], [1])])

end

@testset "TropicalPointConfiguration{min}" begin
  TT = tropical_semiring(min)
  pts1_a = TT[0 -4 -1; 0 -1 0]
  C1_a = tropical_point_configuration(pts1_a)
  @test C1_a isa TropicalPointConfiguration{typeof(min)}

  P1_a = tropical_convex_hull(C1_a)
  @test pm_object(C1_a) == pm_object(P1_a)
end
