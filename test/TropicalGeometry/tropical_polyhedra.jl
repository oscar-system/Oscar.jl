polymake_version = pkgversion(Polymake.polymake_jll)

for minOrMax in (min,max)
__sign = minOrMax == min ? 1 : -1
oo = minOrMax == min ? inf : -inf

@testset "TropicalPolyhedron{$minOrMax}" begin
  TT = tropical_semiring(minOrMax)
  pts1_a = TT[0 -1__sign 0; 0 -4__sign -1__sign]
  pts1_b = [
   TT.(__sign*[0, -1, 0]),
   TT.(__sign*[0, -4, -1]),
   TT.(__sign*[0, -3, 0])
  ]
  P1_a = tropical_convex_hull(pts1_a)

  @test P1_a isa TropicalPolyhedron{typeof(minOrMax)}
  @test ambient_dim(P1_a) == 2
  @test dim(P1_a) == 1
  @test issetequal(vertices(P1_a), eachrow(pts1_a))
  @test n_vertices(P1_a) == 2
  @test issetequal(pseudovertices(P1_a), pts1_b)
  @test n_pseudovertices(P1_a) == 3
  
  @test is_bounded(P1_a) == true

  @test Oscar.pm_object(P1_a) == Oscar.pm_object(P1_a |> tropical_point_configuration)

  cov1 = maximal_covectors(P1_a)
  @test issetequal(cov1, [IncidenceMatrix([[1], [2], [2]]), IncidenceMatrix([[1], [2], [1]])])

  covdec1_a = covector_decomposition(P1_a)
  @test issetequal(covdec1_a |> vertices,
    Vector{QQFieldElem}[
      __sign*[-1, 0],
      __sign*[-4, -1],
      __sign*[-3, 0]
    ])
  @test issetequal(covector_decomposition(P1_a; dehomogenize_by=2) |> vertices,
    Vector{QQFieldElem}[
      __sign*[1, 1],
      __sign*[4, 3],
      __sign*[3, 3]
    ])
  if minOrMax == min
    @test maximal_polyhedra(IncidenceMatrix, covdec1_a) == IncidenceMatrix([[1,3],[2,3]])
  else
    @test maximal_polyhedra(IncidenceMatrix, covdec1_a) == IncidenceMatrix([[1,2],[1,3]])
  end

  ## Adding redundant points gives same tropical polytope with same properties
  P1_b = tropical_convex_hull(pts1_b)
  @test issetequal(P1_a |> vertices, P1_b |> vertices)
  @test dim(P1_b) == 1

  pts2 = TT[0 1__sign 2__sign; oo 0 3__sign; oo oo 0]
  P2 = tropical_convex_hull(pts2)
  @test ambient_dim(P2) == 2
  @test dim(P2) == 2
  @test issetequal(vertices(P2), eachrow(pts2)) 
  @test n_vertices(P2) == 3 
  @test is_bounded(P2) == false 
  @test issetequal(pseudovertices(P2), [TT.(__sign*[0,1,2]),TT.(__sign*[0,-1,2])]) 
  @test n_pseudovertices(P2) == 2 
  @test maximal_covectors(P2) |> length == 1
  @test maximal_covectors(P2) |> first == IncidenceMatrix([[1], [2], [3]]) broken=polymake_version <= v"400.1400.0+0"


  covdec2 = covector_decomposition(P2)
  @test issetequal(covdec2 |> vertices,
    Vector{QQFieldElem}[
      __sign*[-1, 2],
      __sign*[ 1, 2]
    ])
  @test issetequal(covdec2 |> rays .|> Vector{QQFieldElem},
    Vector{QQFieldElem}[
      __sign*[-1, -1],
      __sign*[ 0, -1]
    ])
  @test maximal_polyhedra(IncidenceMatrix, covdec2) == IncidenceMatrix([[1,2,3,4]])

  pts3 = TT[0 1__sign 4__sign; oo 0 2__sign; oo oo 0]
  P3 = tropical_convex_hull(pts3)
  @test ambient_dim(P3) == 2
  @test dim(P3) == 2
  @test issetequal(vertices(P3), eachrow(pts3)) 
  @test n_vertices(P3) == 3 
  @test issetequal(pseudovertices(P3), [TT.(__sign*[0,1,4]),TT.(__sign*[0,1,3])]) 
  @test n_pseudovertices(P3) == 2 
  @test issetequal(maximal_covectors(P3),
    [
      IncidenceMatrix([[1],[1],[2,3]]),
      IncidenceMatrix([[1],[2],[3]])
    ]) broken=polymake_version <= v"400.1400.0+0"

  covdec3 = covector_decomposition(P3)
  @test issetequal(covdec3 |> vertices,
    Vector{QQFieldElem}[
      __sign*[1, 4],
      __sign*[1, 3]
    ])
  @test issetequal(covdec3 |> rays .|> Vector{QQFieldElem},
    Vector{QQFieldElem}[
      __sign*[-1, -1],
      __sign*[ 0, -1]
    ])
  
  if minOrMax == min
    @test maximal_polyhedra(IncidenceMatrix, covdec3) == IncidenceMatrix([[1,2], [2,3,4]])
  else
    @test maximal_polyhedra(IncidenceMatrix, covdec3) == IncidenceMatrix([[1,2], [1,3,4]])
  end

  pts4 = QQ[0 1 0; 0 4 1; 0 3 3; 0 0 2]*__sign
  P4 = tropical_convex_hull(minOrMax,pts4)
  @test ambient_dim(P3) == 2
  @test dim(P3) == 2
  @test issetequal(vertices(P4), eachrow(TT.(pts4)))
  @test n_vertices(P4) == 4
  @test_broken issetequal(pseudovertices(P3), [eachrow(pts3)...,TT.(__sign*[0,1,3])]) 
  @test_broken n_pseudovertices(P3) == 4

  @test issetequal(maximal_covectors(P4),
    [
      IncidenceMatrix([[2],[1,3,4],[2]]),
      IncidenceMatrix([[3],[4],[1,2]]),
      IncidenceMatrix([[3],[1,4],[2]]),
      IncidenceMatrix([[2,3],[4],[1]])
    ])

  @test issetequal(pseudovertices(P4),
    [
     TT.(__sign*[0,1,0]),
     TT.(__sign*[0,4,1]),
     TT.(__sign*[0,3,3]),
     TT.(__sign*[0,0,2]),
     TT.(__sign*[0,0,0]),
     TT.(__sign*[0,2,1]),
     TT.(__sign*[0,3,1]),
     TT.(__sign*[0,3,2]),
     TT.(__sign*[0,1,3]),
     TT.(__sign*[0,0,1])
    ])
end

@testset "TropicalPointConfiguration{$minOrMax}" begin
  TT = tropical_semiring(minOrMax)
  pts1_a = [TT.(__sign*[0, -1, 0]), TT.(__sign*[0, -4, -1])]
  C1_a = tropical_point_configuration(pts1_a)
  @test issetequal(points(C1_a), pts1_a)
  @test n_points(C1_a) == 2
  @test C1_a isa TropicalPointConfiguration{typeof(minOrMax)}

  ## Converting between TropicalPolyhedron and TropicalPointConfiguration should pass around the same polymake object
  P1_a = tropical_convex_hull(C1_a)
  @test Oscar.pm_object(C1_a) == Oscar.pm_object(P1_a)

  cov1 = maximal_covectors(C1_a)
  @test issetequal(cov1,
    [
     IncidenceMatrix([Int[], Int[], [1,2]]),
     IncidenceMatrix([Int[], [2], [1]]),
     IncidenceMatrix([Int[], [2,1], Int[]]),
     IncidenceMatrix([[1], [2], Int[]]),
     IncidenceMatrix([[2,1], Int[], Int[]]),
     IncidenceMatrix([[1], Int[], [2]]),
    ]) 

  ## Making a pseudovertices an input point of the above point configuration
  pts1_b = QQ[0 -1 0; 0 -4 -1; 0 -3 0]*__sign
  C1_b = tropical_point_configuration(minOrMax, pts1_b)
  @test issetequal(points(C1_b), eachrow(TT.(pts1_b)))

  cov1_b = maximal_covectors(C1_b)
  @test issetequal(cov1_b,
    [
     IncidenceMatrix([Int[], Int[], [1,2,3]]),
     IncidenceMatrix([Int[], [2, 3], [1]]),
     IncidenceMatrix([Int[], [1,2,3], Int[]]),
     IncidenceMatrix([[1], [2,3], Int[]]),
     IncidenceMatrix([[1,3], [2], Int[]]),
     IncidenceMatrix([[1,2,3], Int[], Int[]]),
     IncidenceMatrix([[1,3], Int[], [2]])
    ])

  ## Checking the covector decomposition of the torus as `PolyhedralComplex`
  covdec1_a = covector_decomposition(C1_a)
  @test issetequal(covdec1_a |> vertices,
    Vector{QQFieldElem}[
      __sign*[-1, 0],
      __sign*[-4, -1],
      __sign*[-3, 0]
    ])
  @test issetequal(covdec1_a |> rays .|> Vector{QQFieldElem},
    Vector{QQFieldElem}[
      __sign*[1,1],
      __sign*[-1,0],
      __sign*[0,-1]
    ])

  pts2 = TT[0 1__sign 2__sign; oo 0 3__sign; oo oo 0]
  C2 = tropical_point_configuration(pts2)
  @test issetequal(points(C2), eachrow(pts2))
  
  cov2 = maximal_covectors(C2)
  @test issetequal(cov2,
    [
      IncidenceMatrix([[1],[2],[3]]),
      IncidenceMatrix([Int[],[1,2],[3]]),
      IncidenceMatrix([Int[],[2],[1,3]]),
      IncidenceMatrix([[1],Int[],[2,3]]),
      IncidenceMatrix([Int[],Int[],[1,2,3]])
    ]) broken=polymake_version <= v"400.1400.0+0"

  covdec2 = covector_decomposition(C2)
  @test issetequal(covdec2 |> vertices,
    Vector{QQFieldElem}[
      __sign*[-1, 2],
      __sign*[ 1, 2]
    ])
  @test issetequal(covdec2 |> rays .|> Vector{QQFieldElem},
    Vector{QQFieldElem}[
     __sign*[1, 1],
     __sign*[0, -1],
     __sign*[-1, -1],
     __sign*[-1, 0]
    ])

  pts3 = TT[0 1__sign 4__sign; oo 0 2__sign; oo oo 0]
  C3 = tropical_point_configuration(pts3)
  @test issetequal(points(C3), eachrow(pts3))

  @test issetequal(maximal_covectors(C3),
    [
      IncidenceMatrix([Int[],[1],[2,3]]),
      IncidenceMatrix([Int[],Int[],[1,2,3]]),
      IncidenceMatrix([Int[],[1,2],[3]]),
      IncidenceMatrix([[1],Int[],[2,3]]),
      IncidenceMatrix([[1],[2],[3]])
    ]) broken=polymake_version <= v"400.1400.0+0"
end
end
