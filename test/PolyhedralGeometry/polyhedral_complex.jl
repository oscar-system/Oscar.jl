@testset "PolyhedralComplex{$T}" for (f, T) in _prepare_scalar_types()

  I = IncidenceMatrix([[1, 2, 3], [2, 4]])
  P = f.([0 0; 1 0; 0 1; 1 1])
  P2 = f.([0 0 0; 1 0 0; 0 1 0; 1 1 0])
  F = [4]
  L = f.([0 0 1])

  @testset "constructors" begin

    @test polyhedral_complex(f, I, P) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P, F) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P2, F, L) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P2, F, L; non_redundant = true) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P2, nothing, L) isa PolyhedralComplex{T}

  end

  PC = polyhedral_complex(f, I, P)
  PCF = polyhedral_complex(f, I, -P, F)
  PCFL = polyhedral_complex(f, I, P2, F, L)
  PCFLN = polyhedral_complex(f, I, P2, F, L; non_redundant = true)
  PCL = polyhedral_complex(f, I, P2, nothing, L)

  # check non-redundant flag
  @test !Polymake.exists(Oscar.pm_object(PCFLN), "POINTS")
  @test Polymake.exists(Oscar.pm_object(PCFLN), "VERTICES")

  @test Polymake.exists(Oscar.pm_object(PCFL), "POINTS")
  @test !Polymake.exists(Oscar.pm_object(PCFL), "VERTICES")

  @test common_refinement(PC, PCF) isa PolyhedralComplex{T}
  PCR = common_refinement(PC, PCF)

  @test k_skeleton(PC, 1) isa PolyhedralComplex{T}
  PCK = k_skeleton(PC, 1)

  # test constructor with re-arranged arguments
  let PCF2 = polyhedral_complex(f, -P[1:3, :], I[:, 1:3], -P[4:4, :], I[:, 4:4])
    @test vertices(PCF) == vertices(PCF2)
    @test rays(PCF) == rays(PCF2)
    @test IncidenceMatrix(maximal_polyhedra(PCF)) == IncidenceMatrix(maximal_polyhedra(PCF2))
  end

  @testset "core functionality" begin

    @test ambient_dim(PC) == 2
    @test vertices(PC) isa SubObjectIterator{PointVector{T}}
    @test length(vertices(PC)) == 4
    @test point_matrix(vertices(PC)) == matrix(f, P)
    @test vertices(PC) == [[0, 0], [1, 0], [0, 1], [1, 1]]
    @test rays(PCF) isa SubObjectIterator{RayVector{T}}
    @test length(rays(PCF)) == 1
    @test rays(PCF) == [[-1, -1]]
    @test vector_matrix(rays(PCF)) == matrix(f, [-1 -1])
    @test vertices_and_rays(PCFL) isa SubObjectIterator{Union{RayVector{T}, PointVector{T}}}
    @test length(vertices_and_rays(PCFL)) == 4
    @test vertices_and_rays(PCFL) == [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
    @test typeof.(vertices_and_rays(PCFL)) == [PointVector{T}, PointVector{T}, PointVector{T}, RayVector{T}]
    @test vector_matrix(vertices_and_rays(PCFL)) == matrix(f, P2)
    @test maximal_polyhedra(PC) isa SubObjectIterator{Polyhedron{T}}
    @test length(maximal_polyhedra(PC)) == 2
    @test maximal_polyhedra(PC) == convex_hull.([f], [P[1:3, :], P[[2, 4], :]])
    @test number_of_maximal_polyhedra(PC) == 2
    @test is_simplicial(PC)
    @test !is_pure(PCL)
    @test dim(PCL) == 3
    @test polyhedra_of_dim(PC, 1) isa SubObjectIterator{Polyhedron{T}}
    @test length(polyhedra_of_dim(PC, 1)) == 4
    if T == QQFieldElem
      @test polyhedra_of_dim(PC, 1) == convex_hull.(T, [P[[2, 4], :], P[[1, 3], :], P[[1, 2], :], P[[2, 3], :]])
    else
      @test polyhedra_of_dim(PC, 1) == convex_hull.([f], [P[[2, 4], :], P[[1, 3], :], P[[2, 3], :], P[[1, 2], :]])
    end
    @test lineality_space(PCL) isa SubObjectIterator{RayVector{T}}
    @test length(lineality_space(PCL)) == 1
    @test lineality_space(PCL) == [L[:]]
    @test generator_matrix(lineality_space(PCL)) == matrix(QQ, L)

    @test lineality_dim(PCFL) == 1
    @test f_vector(PCL) == [0, 4, 4, 1]
    # Since there is lineality, there are no rays or vertices
    @test n_rays(PCFL) == 0
    @test n_vertices(PCFL) == 0
    @test number_of_polyhedra(PCL) == 9
    @test codim(PCF) == 0
    @test is_embedded(PC)

    mfPCFLN = minimal_faces(PCFLN)
    @test mfPCFLN.base_points == [P2[i, :] for i in 1:3]
    rmlPCFLN = rays_modulo_lineality(PCFLN)
    @test rmlPCFLN.rays_modulo_lineality == [P2[4, :]]
    @test lineality_space(PCFLN) == [L[1, :]]
    @test vertex_indices(maximal_polyhedra(PCFLN)) == I[:, 1:3]
    @test ray_indices(maximal_polyhedra(PCFLN)) == I[:, 4:4]
    @test vertex_and_ray_indices(maximal_polyhedra(PCFLN)) == I
    @test IncidenceMatrix(maximal_polyhedra(PCFLN)) == I
    @test maximal_polyhedra(IncidenceMatrix, PCFLN) == I

    @test polyhedral_complex(maximal_polyhedra(PCF)) isa PolyhedralComplex
    # this should just deepcopy the object
    PCFc = polyhedral_complex(maximal_polyhedra(PCF))
    @test Polymake.list_properties(Oscar.pm_object(PCFc)) ==
            Polymake.list_properties(Oscar.pm_object(PCF))

    c = cube(f, 3, f(-5), f(2))
    # construct from a generic list of polytopes
    @test polyhedral_complex(faces(c, 1)) isa PolyhedralComplex
    @test f_vector(polyhedral_complex(faces(c, 1))) == [8, 12]
    @test f_vector(polyhedral_complex(c)) == [8, 12, 6, 1]

    @test polyhedral_complex(faces(c, 1); non_redundant=true) isa PolyhedralComplex
    @test f_vector(polyhedral_complex(faces(c, 1))) == [8, 12]
    @test f_vector(polyhedral_complex(c)) == [8, 12, 6, 1]

  end

end
