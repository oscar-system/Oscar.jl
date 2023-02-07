@testset "Polyhedral objects with lineality" begin
    @testset "Cone" begin
        C = Cone(fmpq, [0 1], [1 0])
        @test lineality_dim(C) == 1
        @test nrays(C) == 0
        @test length(rays(C)) == 0
        @test length(rays_modulo_lineality(C)) == 1
    end
    
    @testset "Polyhedron" begin
        P = convex_hull(fmpq, [0 0 1], [0 1 0], [1 0 0])
        @test lineality_dim(P) == 1
        @test nrays(P) == 0
        @test nvertices(P) == 0
        @test length(rays(P)) == 0
        @test length(vertices(P)) == 0
        @test length(minimal_faces(P)) == 1
        @test length(rays_modulo_lineality(P)) == 1
    end

    @testset "PolyhedralFan" begin
        P = convex_hull(fmpq, [0 0; 1 0])
        NF = normal_fan(P)
        @test lineality_dim(NF) == 1
        @test nrays(NF) == 0
        @test length(rays(NF)) == 0
        RML = rays_modulo_lineality(NF)
        @test length(RML) == 2
        @test haskey(RML, :lineality_basis)
        @test haskey(RML, :rays_modulo_lineality)
        @test length(RML[:rays_modulo_lineality]) == 2
        @test length(RML[:lineality_basis]) == 1
    end
end
