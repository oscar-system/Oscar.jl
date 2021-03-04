@testset "Group" begin
    
    @testset "Polytope From Group Orbit" begin
        G = symmetric_group(4)
        x = [0,1,2,3]
        M = matrix(ZZ, [permuted(x,g) for g in G])
        P = convex_hull(M)
        @test ambient_dim(P) == 4

        F = facets(P; as = :polyhedra)
        @test nvertices.(F) == [6, 6, 4, 6, 4, 4, 6, 4, 6, 6, 4, 6, 6, 4]

        op = orbit_polytope(x, G)
        @test P == op

    end

    @testset "linear_symmetries" begin
        C0 = cube(2)
        G0 = linear_symmetries(C0)
        @test degree(G0) == 4
    end
end
