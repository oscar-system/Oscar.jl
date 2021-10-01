#TODO: include more examples with nontrivial lineality space

@testset "Polyhedron" begin

    pts = [1 0 0; 0 0 1]'
    @test convex_hull(pts) isa Polyhedron
    Q0 = convex_hull(pts)
    @test convex_hull(pts; non_redundant = true) == Q0
    Q1 = convex_hull(pts, [1 1])
    Q2 = convex_hull(pts, [1 1], [1 1])
    square = cube(2)
    C1 = cube(2, 0, 1)
    Pos = Polyhedron([-1 0 0; 0 -1 0; 0 0 -1], [0,0,0])
    # TODO
    # @test Polyhedron(([-1 0 0; 0 -1 0; 0 0 -1], [0,0,0]); non_redundant = true) == Pos
    L = Polyhedron([-1 0 0; 0 -1 0], [0,0])
    point = convex_hull([0 1 0])
    s = simplex(2)

    @testset "core functionality" begin
        @test nvertices(Q0) == 3
        @test nvertices.(faces(Q0,1)) == [2,2,2]
        @test lattice_points(Q0) isa VectorIterator{PointVector{Polymake.Integer}}
        @test lattice_points(Q0).m == [0 0; 0 1; 1 0]
        @test isfeasible(Q0)
        @test issmooth(Q0)
        @test isnormal(Q0)
        @test isbounded(Q0)
        @test isfulldimensional(Q0)
        @test f_vector(Q0) == [3,3]
        @test intersect(Q0, Q0) == Q0
        @test Q0+Q0 == minkowski_sum(Q0, Q0)
        @test f_vector(Pos) == [1,3,3]
        @test f_vector(L) == [0, 1, 2]
        @test codim(square) == 0
        @test codim(point) == 3
        @test !isfulldimensional(point)
        @test nrays(recession_cone(Pos)) == 3
        @test vertices(point) isa VectorIterator{PointVector{Polymake.Rational}}
        @test vertices(2*point).m == [0 2 0]
        @test vertices([0,1,0] + point).m == [0 2 0]
        @test rays(Pos) isa VectorIterator{RayVector{Polymake.Rational}}
        @test rays(Pos).m == [1 0 0; 0 1 0; 0 0 1]
        @test rays(RayVector, Pos) isa VectorIterator{RayVector{Polymake.Rational}}
        @test lineality_space(L) isa VectorIterator{RayVector{Polymake.Rational}}
        @test lineality_space(L).m == [0 0 1]
        @test faces(square, 1) isa PolyhedronOrConeIterator{Polyhedron}
        @test length(faces(square, 1)) == 4
        @test size(faces(square, 1).lineality) == (0, 3)
        @test isnothing(faces(square, -1))
        @test vertices(minkowski_sum(Q0, square)).m == [2 -1; 2 1; -1 -1; -1 2; 1 2]
        @test facets(Halfspace, Pos) isa HalfspaceIterator{Halfspace}
        @test facets(Pair, Pos) isa HalfspaceIterator{Pair{Polymake.Matrix{Polymake.Rational}, Polymake.Rational}}
        @test facets(Pos) isa HalfspaceIterator{Halfspace}
        @test affine_hull(point) isa HalfspaceIterator{Hyperplane}
        @test affine_hull(point).A == [1 0 0; 0 1 0; 0 0 1] && affine_hull(point).b == [0, 1, 0]
        @test nfacets(square) == 4
        @test lineality_dim(Q0) == 0
        @test nrays(Q1) == 1
        @test lineality_dim(Q2) == 1
    end

    @testset "linear programs" begin
        LP1 = LinearProgram(square,[1,3])
        LP2 = LinearProgram(square,[2,2]; k=3, convention = :min)
        LP3 = LinearProgram(Pos,[1,2,3])

        @test solve_lp(LP1)==(4,[1,1])
        @test solve_lp(LP2)==(-1,[-1,-1])
        @test string(solve_lp(LP3))=="(inf, nothing)"
    end

    @testset "volume" begin
        @test volume(square) == 4
        @test normalized_volume(square) == 8
        @test normalized_volume(s) == 1
    end

    @testset "standard_constructions" begin
        p = upper_bound_theorem(4,8)
        @test p.pm_polytope.F_VECTOR == [8, 28, 40, 20]
        A = archimedean_solid("cuboctahedron")
        @test count(F -> nvertices(F) == 3, faces(A, 2)) == 8
        # due to GLIBCXX issues with the outdated julia-shipped libstdc++
        # we run this only where recent CompilerSupportLibraries are available
        if VERSION >= v"1.6"
            C = catalan_solid("triakis_tetrahedron")
            @test count(F -> nvertices(F) == 3, faces(C, 2)) == 12
        end
        nc = normal_cone(square, 1)
        @test rays(nc).m == [1 0; 0 1]
        @test Polyhedron(facets(A)) == A
    end

end
