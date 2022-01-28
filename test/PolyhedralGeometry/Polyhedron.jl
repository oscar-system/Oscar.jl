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
    L = Polyhedron([-1 0 0; 0 -1 0], [0,0])
    point = convex_hull([0 1 0])
    s = simplex(2)

    @testset "core functionality" begin
        @test nvertices(Q0) == 3
        @test nvertices.(faces(Q0,1)) == [2,2,2]
        @test lattice_points(Q0) isa SubObjectIterator{PointVector{Polymake.Integer}}
        @test point_matrix(lattice_points(Q0)) == matrix(QQ, [0 0; 0 1; 1 0])
        @test matrix(ZZ, lattice_points(Q0)) == matrix(ZZ, [0 0; 0 1; 1 0])
        @test length(lattice_points(Q0)) == 3
        @test lattice_points(Q0)[1] == PointVector{Polymake.Integer}([0, 0])
        @test lattice_points(Q0)[2] == PointVector{Polymake.Integer}([0, 1])
        @test lattice_points(Q0)[3] == PointVector{Polymake.Integer}([1, 0])
        @test interior_lattice_points(square) isa SubObjectIterator{PointVector{Polymake.Integer}}
        @test point_matrix(interior_lattice_points(square)) == matrix(ZZ, [0 0])
        @test matrix(ZZ, interior_lattice_points(square)) == matrix(ZZ,[0 0])
        @test length(interior_lattice_points(square)) == 1
        @test interior_lattice_points(square)[] == PointVector{Polymake.Integer}([0, 0])
        @test boundary_lattice_points(square) isa SubObjectIterator{PointVector{Polymake.Integer}}
        @test point_matrix(boundary_lattice_points(square)) == matrix(ZZ, [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1])
        @test length(boundary_lattice_points(square)) == 8
        @test boundary_lattice_points(square)[1] == PointVector{Polymake.Integer}([-1, -1])
        @test boundary_lattice_points(square)[2] == PointVector{Polymake.Integer}([-1, 0])
        @test boundary_lattice_points(square)[3] == PointVector{Polymake.Integer}([-1, 1])
        @test boundary_lattice_points(square)[4] == PointVector{Polymake.Integer}([0, -1])
        @test boundary_lattice_points(square)[5] == PointVector{Polymake.Integer}([0, 1])
        @test boundary_lattice_points(square)[6] == PointVector{Polymake.Integer}([1, -1])
        @test boundary_lattice_points(square)[7] == PointVector{Polymake.Integer}([1, 0])
        @test boundary_lattice_points(square)[8] == PointVector{Polymake.Integer}([1, 1])
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
        @test vertices(PointVector{Polymake.Rational}, point) isa SubObjectIterator{PointVector{Polymake.Rational}}
        @test vertices(PointVector, point) isa SubObjectIterator{PointVector{Polymake.Rational}}
        @test vertices(point) isa SubObjectIterator{PointVector{Polymake.Rational}}
        @test point_matrix(vertices(2*point)) == matrix(QQ, [0 2 0])
        @test point_matrix(vertices([0,1,0] + point)) == matrix(QQ, [0 2 0])
        @test rays(RayVector{Polymake.Rational}, Pos) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test rays(RayVector, Pos) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test rays(Pos) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test vector_matrix(rays(Pos)) == matrix(QQ, [1 0 0; 0 1 0; 0 0 1])
        @test length(rays(Pos)) == 3
        @test rays(Pos)[1] == RayVector([1, 0, 0])
        @test rays(Pos)[2] == RayVector([0, 1, 0])
        @test rays(Pos)[3] == RayVector([0, 0, 1])
        @test lineality_space(L) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test generator_matrix(lineality_space(L)) == matrix(QQ, [0 0 1])
        @test matrix(ZZ, lineality_space(L)) == matrix(ZZ, [0 0 1])
        @test length(lineality_space(L)) == 1
        @test lineality_space(L)[] == RayVector([0, 0, 1])
        @test faces(square, 1) isa SubObjectIterator{Polyhedron}
        @test length(faces(square, 1)) == 4
        @test faces(square, 1)[1] == convex_hull([-1 -1; -1 1])
        @test faces(square, 1)[2] == convex_hull([1 -1; 1 1])
        @test faces(square, 1)[3] == convex_hull([-1 -1; 1 -1])
        @test faces(square, 1)[4] == convex_hull([-1 1; 1 1])
        @test vertex_indices(faces(square, 1)) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
        @test ray_indices(faces(square, 1)) == IncidenceMatrix(4, 0)
        @test faces(Pos, 1) isa SubObjectIterator{Polyhedron}
        @test length(faces(Pos, 1)) == 3
        @test faces(Pos, 1)[1] == convex_hull([0 0 0], [1 0 0])
        @test faces(Pos, 1)[2] == convex_hull([0 0 0], [0 1 0])
        @test faces(Pos, 1)[3] == convex_hull([0 0 0], [0 0 1])
        @test vertex_indices(faces(Pos, 1)) == IncidenceMatrix([[1], [1], [1]])
        @test ray_indices(faces(Pos, 1)) == IncidenceMatrix([[1], [2], [3]])
        @test isnothing(faces(Q2, 0))
        v = vertices(minkowski_sum(Q0, square))
        @test length(v) == 5
        @test v[1] == PointVector([2, -1])
        @test v[2] == PointVector([2, 1])
        @test v[3] == PointVector([-1, -1])
        @test v[4] == PointVector([-1, 2])
        @test v[5] == PointVector([1, 2])
        @test point_matrix(v) == matrix(QQ, [2 -1; 2 1; -1 -1; -1 2; 1 2])
        for T in [AffineHalfspace, Pair{Polymake.Matrix{Polymake.Rational}, Polymake.Rational}, Polyhedron]
            @test facets(T, Pos) isa SubObjectIterator{T}
            @test length(facets(T, Pos)) == 3
            @test affine_inequality_matrix(facets(T, Pos)) == matrix(QQ, Matrix{fmpq}([0 -1 0 0; 0 0 -1 0; 0 0 0 -1]))
            @test halfspace_matrix_pair(facets(T, Pos)).A == matrix(QQ, Matrix{fmpq}([-1 0 0; 0 -1 0; 0 0 -1])) && halfspace_matrix_pair(facets(T, Pos)).b == [0, 0, 0]
            @test vertex_indices(facets(T, Pos)) == IncidenceMatrix([[1], [1], [1]])
            @test ray_indices(facets(T, Pos)) == IncidenceMatrix([[2, 3], [1, 3], [1, 2]])
            @test facets(T, Pos)[1] == T([-1 0 0], 0)
            @test facets(T, Pos)[2] == T([0 -1 0], 0)
            @test facets(T, Pos)[3] == T([0 0 -1], 0)
        end
        @test facets(Pair, Pos) isa SubObjectIterator{Pair{Polymake.Matrix{Polymake.Rational}, Polymake.Rational}}
        @test facets(Pos) isa SubObjectIterator{AffineHalfspace}
        @test facets(Halfspace, Pos) isa SubObjectIterator{AffineHalfspace}
        @test affine_hull(point) isa SubObjectIterator{AffineHyperplane}
        @test affine_equation_matrix(affine_hull(point)) == matrix(QQ, [0 -1 0 0; 1 0 -1 0; 0 0 0 -1])
        @test Oscar.affine_matrix_for_polymake(affine_hull(point)) == [0 -1 0 0; 1 0 -1 0; 0 0 0 -1]
        @test length(affine_hull(point)) == 3
        @test affine_hull(point)[1] == Hyperplane([1 0 0], 0)
        @test affine_hull(point)[2] == Hyperplane([0 1 0], 1)
        @test affine_hull(point)[3] == Hyperplane([0 0 1], 0)
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
        @test upper_bound_f_vector(4,8) == [8, 28, 40, 20]
        @test upper_bound_g_vector(4,8) == [1, 3, 6]
        @test upper_bound_h_vector(4,8) == [1, 4, 10, 4 ,1]
        A = archimedean_solid("cuboctahedron")
        @test count(F -> nvertices(F) == 3, faces(A, 2)) == 8
        # due to GLIBCXX issues with the outdated julia-shipped libstdc++
        # we run this only where recent CompilerSupportLibraries are available
        if VERSION >= v"1.6"
            C = catalan_solid("triakis_tetrahedron")
            @test count(F -> nvertices(F) == 3, faces(C, 2)) == 12
        end
        nc = normal_cone(square, 1)
        @test vector_matrix(rays(nc)) == matrix(QQ, [1 0; 0 1])
        @test Polyhedron(facets(A)) == A
        b1 = birkhoff(3)
        b2 = birkhoff(3, even = true)
        @test nvertices(pyramid(b1)) + 1 == nvertices(bipyramid(b1))
        @test nvertices(b1) == nvertices(b2) * 2
    end

end
