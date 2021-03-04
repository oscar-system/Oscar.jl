@testset "OscarPolytope" begin

    pts = [1 0 0; 0 0 1]'
    lin = [0 1 0]
    Cone1=positive_hull(pts)
    Q0 = convex_hull(pts)
    Q1 = convex_hull(pts, [1 1])
    Q2 = convex_hull(pts, [1 1], [1 1])
    C0 = cube(2)
    C1 = cube(2, 1, 0)

    @testset "(de)homogenize" begin
        dehomogenize, homogenize = Oscar.dehomogenize, Oscar.homogenize

        m = [1 2 3; 4 5 6]'
        @test dehomogenize(homogenize(m, 0 // 1)) == m
        @test dehomogenize(homogenize(m)) == m
        @test dehomogenize(homogenize(pm.Matrix(m))) == m
        @test dehomogenize(homogenize(pm.Matrix{pm.Integer}(m))) isa pm.Matrix{pm.Integer}
        @test dehomogenize(homogenize(pm.Matrix{pm.Rational}(m))) isa pm.Matrix{pm.Rational}
    end

    @testset "conformance tests" begin
        @test typeof(Cone1) == Cone
        @test typeof(Q0) == Polyhedron
        @test typeof(Q1) == Polyhedron
        @test typeof(Q2) == Polyhedron
        @test typeof(C0) == Polyhedron
        @test typeof(C1) == Polyhedron
        @test typeof(Q0 == Q0) == Bool
        @test typeof(Q0 == Q1) == Bool
        @test typeof(Q0 != Q0) == Bool
        @test typeof(Q0 != Q1) == Bool
        @test Q0 != Q1
        @test C0 != C1
        @test C0 == C0
        @test typeof(dim(Q0)) == Int
        @test typeof(ambient_dim(Q0)) == Int
        @test collect(vertices(Q0)) == collect(vertices(convex_hull(vertices_as_point_matrix(Q0))))
    end

    @testset "convex_hull" begin
        @test size(vertices_as_point_matrix(Q0)) == (3, 2)
        @test size(vertices_as_point_matrix(Q1)) == (3, 2)
        @test size(rays(Q1; as = :point_matrix)) == (1, 2)
        @test size(lineality_space(Q1)) == (0, 2)
        @test size(vertices_as_point_matrix(Q2)) == (2, 2)
        @test size(rays(Q2; as = :point_matrix)) == (0, 2)
        @test size(lineality_space(Q2)) == (1, 2)
        @test dim(Q0) == 2
        @test dim(Q1) == 2
        @test dim(Q2) == 2
        @test ambient_dim(Q0) == 2
        @test ambient_dim(Q1) == 2
        @test ambient_dim(Q2) == 2
    end

    @testset "standard constructions" begin
        @test size(vertices_as_point_matrix(C0)) == (4,2)
        @test C0 == convex_hull(vertices_as_point_matrix(C0))
        @test isbounded(C0)
        @test issmooth(C0)
        @test isnormal(C0)
        @test isfeasible(C0)
        @test isfulldimensional(C0)
        @test minkowski_sum(C0,C0) == cube(2,2,-2)
        @test minkowski_sum(C0,C0; algorithm=:fukuda) == cube(2,2,-2)
        @test intersect(C0,C0) == C0
    end

    @testset "newton_polytope" begin
        Qx, x = Oscar.PolynomialRing(Oscar.QQ, :x => 1:2)
        f = sum([x; 1])^2 + x[1]^4 * x[2] * 3
        newt = newton_polytope(f)
        @test dim(newt) == 2
        @test vertices_as_point_matrix(newt) == [4 1; 2 0; 0 2; 0 0]
    end

    @testset "Construct from fmpq" begin
        A = zeros(Oscar.QQ, 3, 2)
        A[1, 1] = 1
        A[3, 2] = 4
        @test vertices_as_point_matrix(convex_hull(A)) == [1 0; 0 0; 0 4]

        lhs, rhs = facets_as_halfspace_matrix_pair(Polyhedron(A, [1, 2, -3]))
        @test lhs == [1 0; 0 4; 0 0]
        @test rhs == [1, -3, 1]
    end

    @testset "Construct from fmpz" begin
        A = zeros(Oscar.ZZ, 3, 2)
        A[1, 1] = 1
        A[3, 2] = 4
        @test vertices_as_point_matrix(convex_hull(A)) == [1 0; 0 0; 0 4]

        lhs, rhs = facets_as_halfspace_matrix_pair(Polyhedron(A, [1, 2, -3]))
        @test lhs == [1 0; 0 4; 0 0]
        @test rhs == [1, -3, 1]
    end

end # of @testset "OscarPolytope"

