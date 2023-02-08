@testset "Polyhedral objects with lineality" begin
    @testset "Cone" begin
        C = Cone(fmpq, [0 1], [1 0])
        @test lineality_dim(C) == 1
        @test nrays(C) == 0
        @test length(rays(C)) == 0
        
        RML = rays_modulo_lineality(C)
        @test RML isa NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{fmpq}}, SubObjectIterator{RayVector{fmpq}}}}
        @test length(RML) == 2
        @test haskey(RML, :lineality_basis)
        @test haskey(RML, :rays_modulo_lineality)
        @test length(RML[:rays_modulo_lineality]) == 1
        @test length(RML[:lineality_basis]) == 1
    end
    
    @testset "Polyhedron" begin
        P = convex_hull(fmpq, [0 0 1], [0 1 0], [1 0 0])
        @test lineality_dim(P) == 1
        @test nrays(P) == 0
        @test nvertices(P) == 0
        @test length(rays(P)) == 0
        @test length(vertices(P)) == 0
        
        MFP = minimal_faces(P)
        @test MFP isa NamedTuple{(:base_points, :lineality_basis), Tuple{SubObjectIterator{PointVector{fmpq}}, SubObjectIterator{RayVector{fmpq}}}}
        @test length(MFP) == 2
        @test haskey(MFP, :lineality_basis)
        @test haskey(MFP, :base_points)
        @test length(MFP[:base_points]) == 1
        @test length(MFP[:lineality_basis]) == 1
        
        RML = rays_modulo_lineality(P)
        @test RML isa NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{fmpq}}, SubObjectIterator{RayVector{fmpq}}}}
        @test length(RML) == 2
        @test haskey(RML, :lineality_basis)
        @test haskey(RML, :rays_modulo_lineality)
        @test length(RML[:rays_modulo_lineality]) == 1
        @test length(RML[:lineality_basis]) == 1
    end

    @testset "PolyhedralFan" begin
        P = convex_hull(fmpq, [0 0; 1 0])
        NF = normal_fan(P)
        @test lineality_dim(NF) == 1
        @test nrays(NF) == 0
        @test length(rays(NF)) == 0
        
        RML = rays_modulo_lineality(NF)
        @test RML isa NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{fmpq}}, SubObjectIterator{RayVector{fmpq}}}}
        @test length(RML) == 2
        @test haskey(RML, :lineality_basis)
        @test haskey(RML, :rays_modulo_lineality)
        @test length(RML[:rays_modulo_lineality]) == 2
        @test length(RML[:lineality_basis]) == 1
    end

    @testset "PolyhedralComplex" begin
        VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0]
        IM = IncidenceMatrix([[1,2,3],[1,3,4]])
        far_vertices = [2,3,4]
        L = [0 0 1]
        PC = PolyhedralComplex(IM, VR, far_vertices, L)
        @test length(vertices(PC)) == 0
        @test length(rays(PC)) == 0
        
        MFP = minimal_faces(PC)
        @test MFP isa NamedTuple{(:base_points, :lineality_basis), Tuple{SubObjectIterator{PointVector{fmpq}}, SubObjectIterator{RayVector{fmpq}}}}
        @test length(MFP) == 2
        @test haskey(MFP, :lineality_basis)
        @test haskey(MFP, :base_points)
        @test length(MFP[:base_points]) == 1
        @test length(MFP[:lineality_basis]) == 1
        
        RML = rays_modulo_lineality(PC)
        @test RML isa NamedTuple{(:rays_modulo_lineality, :lineality_basis), Tuple{SubObjectIterator{RayVector{fmpq}}, SubObjectIterator{RayVector{fmpq}}}}
        @test length(RML) == 2
        @test haskey(RML, :lineality_basis)
        @test haskey(RML, :rays_modulo_lineality)
        @test length(RML[:rays_modulo_lineality]) == 3
        @test length(RML[:lineality_basis]) == 1
    end
end
