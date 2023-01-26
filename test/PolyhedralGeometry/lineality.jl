@testset "Polyhedral objects with lineality" begin
    @testset "Cone" begin
        C = Cone(fmpq, [0 1], [1 0])
        @test lineality_dim(C) == 1
        @test nrays(C) == 0
    end
    
    @testset "Polyhedron" begin
        P = convex_hull(fmpq, [0 0 1], [0 1 0], [1 0 0])
        @test lineality_dim(P) == 1
        @test nrays(P) == 0
        @test nvertices(P) == 0
    end
end
