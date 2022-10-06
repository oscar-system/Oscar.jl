@testset "Group" begin

    @testset "Polytope From Group Orbit" begin
        G = symmetric_group(4)
        x = [0,1,2,3]
        M = matrix(ZZ, [permuted(x,g) for g in G])
        P = convex_hull(M)
        @test ambient_dim(P) == 4

        F = facets(Polyhedron, P)
        #@test nvertices.(F) == [6, 6, 4, 6, 4, 4, 6, 4, 6, 6, 4, 6, 6, 4]
        #Since different convex hull algorithms will result in differnt vertex orders, this test may fail.
        #Best way to handle this is to use "prefer" option in polymake, which is not available in OSCAR yet
        #We avoid using "sort"(expensive) or "countmap" from StatsBase(unneccesary dependency), and just check the lenghts for now.
        @test length(nvertices.(F)) == 14

        op = orbit_polytope(x, G)
        @test P == op

    end

    @testset "linear_symmetries" begin
        C0 = cube(2)
        G0 = linear_symmetries(C0)
        @test degree(G0) == 4
        P = Polyhedron([-1 0 0; 0 -1 0; 0 0 -1],[0,0,0])
        @test_throws ArgumentError linear_symmetries(P)
    end

    @testset "vf_group" begin
        C = cube(3)
        A = vf_group(C)
        @test degree(A) == 6
        P = Polyhedron([-1 0 0; 0 -1 0; 0 0 -1],[0,0,0])
        @test_throws ArgumentError vf_group(P)
    end

    @testset "combinatorial_symmetries" begin
        C = cube(3)
        A = combinatorial_symmetries(C)
        @test degree(A) == 8
        P = Polyhedron([-1 0 0; 0 -1 0; 0 0 -1],[0,0,0])
        @test_throws ArgumentError combinatorial_symmetries(P)
    end
end
