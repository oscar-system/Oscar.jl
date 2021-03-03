@testset "Polyhedron" begin
    
    pts = [1 0 0; 0 0 1]'
    Q0 = convex_hull(pts)
    Q1 = convex_hull(pts, [1 1])
    Q2 = convex_hull(pts, [1 1], [1 1])
    C0 = cube(2)
    C1 = cube(2, 1, 0)
    Pos=Polyhedron([-1 0 0; 0 -1 0; 0 0 -1],[0,0,0])
    point = convex_hull([0 1 0])
    s = simplex(2)

    @testset "core functionality" begin
        @test nvertices(Q0) == 3
        @test nvertices.(faces(Q0,1)) == [2,2,2]
        @test length.(lattice_points(Q0)) == [2;2;2]
        @test isfeasible(Q0)
        @test issmooth(Q0)
        @test isnormal(Q0)
        @test isbounded(Q0)
        @test isfulldimensional(Q0)
        @test f_vector(Q0) == [3,3]
        @test intersect(Q0, Q0) == Q0
        @test Q0+Q0 == minkowski_sum(Q0, Q0)
        @test f_vector(Pos) == [1,3,3]
        @test codim(C0) == 0
        @test codim(point) == 3
        @test !isfulldimensional(point)
        @test nrays(recession_cone(Pos)) == 3
        @test vertices_as_point_matrix(2*point) == [0 2 0]
        @test vertices_as_point_matrix([0,1,0] + point) == [0 2 0]
        @test length(collect(rays(Pos))) == 3
    end
    
    @testset "linear programs" begin
        LP1 = LinearProgram(C0,[1,3])
        LP2 = LinearProgram(C0,[2,2],3; convention = :min)
        LP3 = LinearProgram(Pos,[1,2,3])

        @test solve_lp(LP1)==(4,[1,1])
        @test solve_lp(LP2)==(-1,[-1,-1])
        @test string(solve_lp(LP3))=="(inf, nothing)"
    end

    @testset "volume" begin
        @test volume(C0) == 4
        @test normalized_volume(C0) == 8
        @test normalized_volume(s) == 1
    end

end
