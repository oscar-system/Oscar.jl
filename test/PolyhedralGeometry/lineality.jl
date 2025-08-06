@testset "Polyhedral objects with lineality" begin
  @testset "Cone" begin
    C = positive_hull(QQFieldElem, [0 1], [1 0])
    @test lineality_dim(C) == 1
    @test n_rays(C) == 0
    @test length(rays(C)) == 0
    @test size(vector_matrix(rays(C))) == (0, 2)

    RML = rays_modulo_lineality(C)
    @test RML isa NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{
        SubObjectIterator{RayVector{QQFieldElem}},SubObjectIterator{RayVector{QQFieldElem}}
      },
    }
    @test length(RML) == 2
    @test haskey(RML, :lineality_basis)
    @test haskey(RML, :rays_modulo_lineality)
    @test length(RML[:rays_modulo_lineality]) == 1
    @test length(RML[:lineality_basis]) == 1
  end

  @testset "Polyhedron" begin
    P = convex_hull(QQFieldElem, [0 0 1], [0 1 0], [1 0 0])
    @test lineality_dim(P) == 1
    @test n_rays(P) == 0
    @test n_vertices(P) == 0
    @test length(rays(P)) == 0
    @test size(vector_matrix(rays(P))) == (0, 3)
    @test length(vertices(P)) == 0
    @test size(point_matrix(vertices(P))) == (0, 3)

    MFP = minimal_faces(P)
    @test MFP isa NamedTuple{
      (:base_points, :lineality_basis),
      Tuple{
        SubObjectIterator{PointVector{QQFieldElem}},
        SubObjectIterator{RayVector{QQFieldElem}},
      },
    }
    @test length(MFP) == 2
    @test haskey(MFP, :lineality_basis)
    @test haskey(MFP, :base_points)
    @test length(MFP[:base_points]) == 1
    @test length(MFP[:lineality_basis]) == 1

    RML = rays_modulo_lineality(P)
    @test RML isa NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{
        SubObjectIterator{RayVector{QQFieldElem}},SubObjectIterator{RayVector{QQFieldElem}}
      },
    }
    @test length(RML) == 2
    @test haskey(RML, :lineality_basis)
    @test haskey(RML, :rays_modulo_lineality)
    @test length(RML[:rays_modulo_lineality]) == 1
    @test length(RML[:lineality_basis]) == 1
  end

  @testset "PolyhedralFan" begin
    P = convex_hull(QQFieldElem, [0 0; 1 0])
    NF = normal_fan(P)
    @test lineality_dim(NF) == 1
    @test n_rays(NF) == 0
    @test length(rays(NF)) == 0
    @test size(vector_matrix(rays(NF))) == (0, 2)

    RML = rays_modulo_lineality(NF)
    @test RML isa NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{
        SubObjectIterator{RayVector{QQFieldElem}},SubObjectIterator{RayVector{QQFieldElem}}
      },
    }
    @test length(RML) == 2
    @test haskey(RML, :lineality_basis)
    @test haskey(RML, :rays_modulo_lineality)
    @test length(RML[:rays_modulo_lineality]) == 2
    @test length(RML[:lineality_basis]) == 1
  end

  @testset "PolyhedralComplex" begin
    VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0]
    IM = incidence_matrix([[1, 2, 3], [1, 3, 4]])
    far_vertices = [2, 3, 4]
    L = [0 0 1]
    PC = polyhedral_complex(IM, VR, far_vertices, L)
    @test length(vertices(PC)) == 0
    @test size(point_matrix(vertices(PC))) == (0, 3)
    @test length(rays(PC)) == 0
    @test size(vector_matrix(rays(PC))) == (0, 3)

    MFP = minimal_faces(PC)
    @test MFP isa NamedTuple{
      (:base_points, :lineality_basis),
      Tuple{
        SubObjectIterator{PointVector{QQFieldElem}},
        SubObjectIterator{RayVector{QQFieldElem}},
      },
    }
    @test length(MFP) == 2
    @test haskey(MFP, :lineality_basis)
    @test haskey(MFP, :base_points)
    @test length(MFP[:base_points]) == 1
    @test length(MFP[:lineality_basis]) == 1

    RML = rays_modulo_lineality(PC)
    @test RML isa NamedTuple{
      (:rays_modulo_lineality, :lineality_basis),
      Tuple{
        SubObjectIterator{RayVector{QQFieldElem}},SubObjectIterator{RayVector{QQFieldElem}}
      },
    }
    @test length(RML) == 2
    @test haskey(RML, :lineality_basis)
    @test haskey(RML, :rays_modulo_lineality)
    @test length(RML[:rays_modulo_lineality]) == 3
    @test length(RML[:lineality_basis]) == 1
  end
end
