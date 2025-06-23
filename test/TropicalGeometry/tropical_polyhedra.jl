@testset "TropicalPolyhedron{min}" begin
  TT = tropical_semiring(min)
  pts1_a = TT[0 -1 0; 0 -4 -1]
  P1_a = tropical_convex_hull(pts1_a)
  @test P1_a isa TropicalPolyhedron{typeof(min)}
  @test dim(P1_a) == 1

  ## Adding redundant points gives same tropical polytope with same properties
  pts1_b = TT[0 -1 0; 0 -4 -1; 0 -3 0]
  P1_b = tropical_convex_hull(pts1_b)
  @test issetequal(P1_a |> vertices, P1_b |> vertices)
  @test dim(P1_b) == 1

  cov1 = maximal_covectors(P1_a)
  @test issetequal(cov1, [IncidenceMatrix([[1], [2], [2]]), IncidenceMatrix([[1], [2], [1]])])

  pts2 = TT[0 1 2; inf 0 3; inf inf 0]
  P2 = tropical_convex_hull(pts2)
  @test dim(P2) == 2
  @test maximal_covectors(P2) |> length == 1
  @test maximal_covectors(P2) |> first == IncidenceMatrix([[1], [2], [3]])
end

@testset "TropicalPointConfiguration{min}" begin
  TT = tropical_semiring(min)
  pts1_a = TT[0 -1 0; 0 -4 -1]
  C1_a = tropical_point_configuration(pts1_a)
  @test C1_a isa TropicalPointConfiguration{typeof(min)}

  P1_a = tropical_convex_hull(C1_a)
  @test pm_object(C1_a) == pm_object(P1_a)

  cov1 = maximal_covectors(C1_a)
  @test issetequal(cov1,
    [
     IncidenceMatrix([Int[], Int[], [1,2]]),
     IncidenceMatrix([Int[], [1], [2]]),
     IncidenceMatrix([Int[], [1,2], Int[]]),
     IncidenceMatrix([[2], [1], Int[]]),
     IncidenceMatrix([[1,2], Int[], Int[]]),
     IncidenceMatrix([[2], Int[], [1]]),
    ])

  pts1_b = TT[0 -1 0; 0 -4 -1; 0 -3 0]
  C1_b = tropical_convex_hull(pts1_b)
  cov1_b = maximal_covectors(C1_b)
  @test issetequal(cov1_b,
    [
     IncidenceMatrix([Int[], Int[], [1,2,3]]),
     IncidenceMatrix([Int[], [1, 3], [2]]),
     IncidenceMatrix([Int[], [1,2,3], Int[]]),
     IncidenceMatrix([[2], [1,3], Int[]]),
     IncidenceMatrix([[2,3], [1], Int[]]),
     IncidenceMatrix([[1,2,3], Int[], Int[]]),
     IncidenceMatrix([[2,3], Int[], [1]])
    ])

  covdec1_a = covector_decomposition(C1_a)
  @test issetequal(covdec1_a |> vertices,
    [
      [-1, 0],
      [-4, -1],
      [-3, 0]
    ])
  @test issetequal(covdec1_a |> rays,
    [
      [1,1],
      [-1,0],
      [0,-1]
    ])

  pts2 = TT[0 1 2; inf 0 3; inf inf 0]
  C2 = tropical_point_configuration(pts2)
  cov2 = maximal_covectors(C2)
  @test issetequal(cov2,
    [
      IncidenceMatrix([[1],[2],[3]]),
      IncidenceMatrix([Int[],[12],[3]]),
      IncidenceMatrix([Int[],[2],[13]]),
      IncidenceMatrix([[1],Int[],[23]]),
      IncidenceMatrix([Int[],Int[],[123]])
    ])
end
