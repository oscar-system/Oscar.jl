@testset "PolyhedralComplex{$T}" for (f, T) in _prepare_scalar_types()
  I = incidence_matrix([[1, 2, 3], [2, 4]])
  P = f.([0 0; 1 0; 0 1; 1 1])
  P2 = f.([0 0 0; 1 0 0; 0 1 0; 1 1 0])
  F = [4]
  L = f.([0 0 1])

  @testset "constructors" begin
    @test polyhedral_complex(f, I, P) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P, F) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P2, F, L) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P2, F, L; non_redundant=true) isa PolyhedralComplex{T}
    @test polyhedral_complex(f, I, P2, nothing, L) isa PolyhedralComplex{T}
  end

  PC = polyhedral_complex(f, I, P)
  PCF = polyhedral_complex(f, I, -P, F)
  PCFL = polyhedral_complex(f, I, P2, F, L)
  PCFLN = polyhedral_complex(f, I, P2, F, L; non_redundant=true)
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
    @test issetequal(vertices(PCF), vertices(PCF2))
    @test issetequal(rays(PCF), rays(PCF2))
    @test issetequal(maximal_polyhedra(PCF), maximal_polyhedra(PCF2))
  end

  @testset "core functionality" begin
    @test ambient_dim(PC) == 2
    @test vertices(PC) isa SubObjectIterator{PointVector{T}}
    @test length(vertices(PC)) == 4
    @test issetequal(vertices(PC), point_vector.(Ref(f), [[0, 0], [1, 0], [0, 1], [1, 1]]))
    @test point_matrix(vertices(PC)) == _oscar_matrix_from_property(f, vertices(PC))
    @test rays(PCF) isa SubObjectIterator{RayVector{T}}
    @test length(rays(PCF)) == 1
    @test rays(PCF) == [[-1, -1]]
    @test vector_matrix(rays(PCF)) == matrix(f, [-1 -1])
    @test vertices_and_rays(PCFL) isa SubObjectIterator{Union{RayVector{T},PointVector{T}}}
    @test length(vertices_and_rays(PCFL)) == 4
    @test issetequal(
      vertices_and_rays(PCFL),
      [
        point_vector(f, [0, 0, 0]),
        point_vector(f, [1, 0, 0]),
        point_vector(f, [0, 1, 0]),
        ray_vector(f, [1, 1, 0]),
      ],
    )
    @test vector_matrix(vertices_and_rays(PCFL)) ==
      _oscar_matrix_from_property(f, vertices_and_rays(PCFL))
    @test maximal_polyhedra(PC) isa SubObjectIterator{Polyhedron{T}}
    @test length(maximal_polyhedra(PC)) == 2
    @test issetequal(maximal_polyhedra(PC), convex_hull.([f], [P[1:3, :], P[[2, 4], :]]))
    @test number_of_maximal_polyhedra(PC) == 2
    @test is_simplicial(PC)
    @test !is_pure(PCL)
    @test dim(PCL) == 3
    @test polyhedra_of_dim(PC, 1) isa SubObjectIterator{Polyhedron{T}}
    @test length(polyhedra_of_dim(PC, 1)) == 4
    @test issetequal(polyhedra_of_dim(PC, 1),
      convex_hull.(Ref(f), [P[[2, 4], :], P[[1, 3], :], P[[2, 3], :], P[[1, 2], :]]))
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
    @test incidence_matrix(maximal_polyhedra(PCFLN)) == I
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

  @testset "Fan conversion" begin
    F1 = normal_fan(cube(f, 2))
    F2 = normal_fan(convex_hull(f, [0 0; 1 0]))
    IM = incidence_matrix([[1, 2], [2, 3], [4]])
    R = [0 1 0; 0 0 1; 0 -1 0; 0 -1 -1]
    F3 = polyhedral_fan(f, IM, R, [1 0 0])
    for F in [F1, F2, F3]
      PC = polyhedral_complex(F)
      @test dim(F) == dim(PC)
      @test ambient_dim(F) == ambient_dim(PC)
      @test lineality_dim(F) == lineality_dim(PC)
      @test issetequal(rays(F), rays(PC))
      @test n_maximal_cones(F) == n_maximal_polyhedra(PC)
      vrep = PolyhedralComplex{T}(
        Polymake.fan.PolyhedralComplex(;
          INPUT_RAYS=Oscar.pm_object(PC).RAYS,
          INPUT_CONES=Oscar.pm_object(PC).MAXIMAL_CONES,
          INPUT_LINEALITY=Oscar.pm_object(PC).LINEALITY_SPACE,
        ),
        f,
      )
      hrep = PolyhedralComplex{T}(
        Polymake.fan.PolyhedralComplex(;
          FACET_NORMALS=Oscar.pm_object(PC).FACET_NORMALS,
          MAXIMAL_CONES_FACETS=Oscar.pm_object(PC).MAXIMAL_CONES_FACETS,
          LINEAR_SPAN_NORMALS=Oscar.pm_object(PC).LINEAR_SPAN_NORMALS,
          MAXIMAL_CONES_LINEAR_SPAN_NORMALS=Oscar.pm_object(
            PC
          ).MAXIMAL_CONES_LINEAR_SPAN_NORMALS,
        ),
        f,
      )
      @test f_vector(vrep) == f_vector(hrep)
      @test n_maximal_polyhedra(vrep) == n_maximal_polyhedra(hrep)
    end
  end
end
