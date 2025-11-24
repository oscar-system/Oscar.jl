@testset "Conformance tests" begin
  pm = Polymake

  pts = [1 0; 0 0; 0 1]
  lin = [0 1 0]
  Cone1 = positive_hull(pts)
  Q0 = convex_hull(pts)
  Q1 = convex_hull(pts, [1 1])
  Q2 = convex_hull(pts, [1 1], [1 1])
  C0 = cube(2)
  C1 = cube(2, 0, 1)

  @testset "(de)homogenize" begin
    dehomogenize, homogenize = Oscar.dehomogenize, Oscar.homogenize

    m = [1 2; 3 4; 5 6]
    @test dehomogenize(homogenize(QQ, m, 0//1)) == m
    @test dehomogenize(homogenize(QQ, m)) == m
    @test dehomogenize(homogenize(QQ, pm.Matrix(m))) == m
    @test dehomogenize(homogenize(ZZ, pm.Matrix{pm.Integer}(m))) isa pm.Matrix{pm.Integer}
    @test dehomogenize(homogenize(QQ, pm.Matrix{pm.Rational}(m))) isa pm.Matrix{pm.Rational}
  end

  @testset "conformance tests" begin
    @test typeof(Cone1) == Cone{QQFieldElem}
    @test typeof(Q0) == Polyhedron{QQFieldElem}
    @test typeof(Q1) == Polyhedron{QQFieldElem}
    @test typeof(Q2) == Polyhedron{QQFieldElem}
    @test typeof(C0) == Polyhedron{QQFieldElem}
    @test typeof(C1) == Polyhedron{QQFieldElem}
    @test typeof(Q0 == Q0) == Bool
    @test typeof(Q0 == Q1) == Bool
    @test typeof(Q0 != Q0) == Bool
    @test typeof(Q0 != Q1) == Bool
    @test Q0 != Q1
    @test C0 != C1
    @test C0 == C0
    @test typeof(dim(Q0)) == Int
    @test typeof(ambient_dim(Q0)) == Int
    @test collect(vertices(Q0)) == collect(vertices(convex_hull(vertices(Q0))))
  end

  @testset "convex_hull" begin
    @test size(point_matrix(vertices(Q0))) == (3, 2)
    @test size(point_matrix(vertices(Q1))) == (3, 2)
    @test size(vector_matrix(rays(Q1))) == (1, 2)
    @test size(generator_matrix(lineality_space(Q1))) == (0, 2)
    @test size(point_matrix(vertices(Q2))) == (0, 2)
    @test size(vector_matrix(rays(Q2))) == (0, 2)
    @test size(generator_matrix(lineality_space(Q2))) == (1, 2)
    @test dim(Q0) == 2
    @test dim(Q1) == 2
    @test dim(Q2) == 2
    @test ambient_dim(Q0) == 2
    @test ambient_dim(Q1) == 2
    @test ambient_dim(Q2) == 2
  end

  @testset "standard constructions" begin
    @test size(point_matrix(vertices(C0))) == (4, 2)
    @test C0 == convex_hull(vertices(C0))
    @test is_bounded(C0)
    @test is_smooth(C0)
    @test is_normal(C0)
    @test is_feasible(C0)
    @test is_fulldimensional(C0)
    @test minkowski_sum(C0, C0) == cube(2, -2, 2)
    @test minkowski_sum(C0, C0; algorithm=:fukuda) == cube(2, -2, 2)
    @test intersect(C0, C0) == C0
    Cs = [C0, C0, C0]
    @test intersect(Cs) isa Polyhedron{QQFieldElem}
    @test intersect(Cs...) isa Polyhedron{QQFieldElem}
    @test intersect(Cs) == C0
    @test intersect(Cs...) == C0
    @test convex_hull(Cs) isa Polyhedron{QQFieldElem}
    @test convex_hull(Cs...) isa Polyhedron{QQFieldElem}
    @test convex_hull(Cs) == C0
    @test convex_hull(Cs...) == C0
  end

  @testset "newton_polytope" begin
    Qx, x = Oscar.polynomial_ring(Oscar.QQ, :x => 1:2)
    f = sum([x; 1])^2 + x[1]^4 * x[2] * 3
    newt = newton_polytope(f)
    @test dim(newt) == 2
    @test issetequal(
      vertices(newt), point_vector.(Ref(QQ), [[4, 1], [2, 0], [0, 2], [0, 0]])
    )
  end

  @testset "Construct from QQFieldElem" begin
    A = zeros(Oscar.QQ, 3, 2)
    A[1, 1] = 1
    A[3, 2] = 4
    @test issetequal(
      vertices(convex_hull(A)), point_vector.(Ref(QQ), [[1, 0], [0, 0], [0, 4]])
    )

    @test issetequal(
      facets(polyhedron(A, [1, 2, -3])),
      [affine_halfspace(QQ, [1, 0], 1), affine_halfspace(QQ, [0, 4], -3)],
    )
  end

  @testset "Construct from ZZRingElem" begin
    A = zeros(Oscar.ZZ, 3, 2)
    A[1, 1] = 1
    A[3, 2] = 4
    @test issetequal(
      vertices(convex_hull(A)), point_vector.(Ref(QQ), [[1, 0], [0, 0], [0, 4]])
    )

    @test issetequal(
      facets(polyhedron(A, [1, 2, -3])),
      [affine_halfspace(QQ, [1, 0], 1), affine_halfspace(QQ, [0, 4], -3)],
    )
  end

  @testset "SubObjectIterator/Matrix compatibility" begin
    Pos_poly = convex_hull([0 0 0], [1 0 0; 0 1 0; 0 0 1])
    Pos_cone = positive_hull([1 0 0; 0 1 0; 0 0 1])

    @test Oscar._ambient_dim([1, 2, 3, 4, 5]) == 5
    @test Oscar._ambient_dim([1 2 3 4; 5 6 7 8]) == 4
    @test Oscar._ambient_dim([[1, 2, 3], [4, 5, 6], [7, 8, 9], [0, 1, 2]]) == 3
    @test Oscar._ambient_dim(vertices(Pos_poly)) == 3
    @test Oscar._ambient_dim(collect(vertices(Pos_poly))) == 3

    # test for correctness of input interpretation for different types

    @test convex_hull(vertices(Pos_poly), rays(Pos_poly)) == Pos_poly
    @test positive_hull(rays(Pos_cone)) == Pos_cone

    @test convex_hull([-1 -1 -1], rays(Pos_cone)) == Pos_poly + [-1, -1, -1]
    @test positive_hull(rays(Pos_poly)) == Pos_cone

    @test convex_hull([[0, 0, 0]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == Pos_poly
    @test positive_hull([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == Pos_cone

    @test convex_hull([0 0 0], matrix(ZZ, [1 0 0; 0 1 0; 0 0 1])) == Pos_poly

    @test convex_hull(collect(vertices(Pos_poly)), collect(rays(Pos_poly))) == Pos_poly
    @test positive_hull(collect(rays(Pos_poly))) == Pos_cone

    @test convex_hull([0, 0, 0], rays(Pos_poly)) == Pos_poly
    @test rays(positive_hull([1, 0, 0]))[] == [1, 0, 0]

    @test polyhedron(facets(Pos_poly)) == Pos_poly
    @test polyhedron(facets(Pos_cone)) == Pos_poly

    @test cone_from_inequalities(facets(Pos_poly)) == Pos_cone
    @test cone_from_inequalities(facets(Pos_cone)) == Pos_cone

    @test polyhedron(collect(facets(Pos_poly))) == Pos_poly
    @test polyhedron(collect(facets(Pos_cone))) == Pos_poly

    @test cone_from_inequalities(collect(facets(Pos_poly))) == Pos_cone
    @test cone_from_inequalities(collect(facets(Pos_cone))) == Pos_cone

    # testing correct dispatch and tuple processing for Polyhedron
    @test polyhedron([-1 0 0; 0 -1 0; 0 0 -1], [0, 0, 0]) == Pos_poly
    @test polyhedron([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], QQFieldElem[0, 0, 0]) == Pos_poly
    @test polyhedron(matrix(ZZ, [-1 0 0; 0 -1 0; 0 0 -1]), [0, 0, 0]) == Pos_poly
    @test polyhedron(matrix(QQ, [-1 0 0; 0 -1 0; 0 0 -1]), [0, 0, 0]) == Pos_poly

    # testing different input types
    @test polyhedron([-1 0 0; 0 -1 0; 0 0 -1], Float64[0, 0, 0]) isa Polyhedron{Float64}
    @test polyhedron(Float64[-1 0 0; 0 -1 0; 0 0 -1], [0, 0, 0]) isa Polyhedron{Float64}

    let y = convex_hull([0, 0, 0], [1, 0, 0], [[0, 1, 0], [0, 0, 1]])
      @test polyhedron([-1 0 0], [0]) == y
      @test polyhedron([-1 0 0], 0) == y
      @test polyhedron([-[1, 0, 0]], QQFieldElem[0]) == y
      @test polyhedron([-1, 0, 0], QQFieldElem[0]) == y
      @test polyhedron([[-1, 0, 0]], QQFieldElem(0)) == y
      @test polyhedron([-1, 0, 0], QQFieldElem(0)) == y
      @test polyhedron(matrix(ZZ, [-1 0 0]), [0]) == y
      @test polyhedron(matrix(QQ, [-1 0 0]), [0]) == y
    end

    let x = positive_hull([1 0 0; 0 1 0]), y = convex_hull([0 0 0], [1 0 0; 0 1 0])
      @test polyhedron(facets(y), affine_hull(y)) == y
      @test polyhedron(facets(y), linear_span(x)) == y

      @test cone_from_inequalities(facets(y), affine_hull(y)) == x
      @test cone_from_inequalities(facets(x), linear_span(x)) == x

      @test polyhedron(facets(y), collect(affine_hull(y))) == y
      @test polyhedron(facets(x), collect(linear_span(x))) == y

      @test cone_from_inequalities(facets(x), collect(affine_hull(y))) == x
      @test cone_from_inequalities(facets(y), collect(linear_span(x))) == x

      @test cone_from_inequalities([[-1, 0, 0], [0, -1, 0]], [0 0 1]) == x
      @test cone_from_inequalities([-1 0 0; 0 -1 0], [[0, 0, 1]]) == x
    end

    # Here the content of the SubObjectIterator does not fit the idea of the
    # methods; we want ArgumentErrors to be thrown
    @test_throws ArgumentError convex_hull(facets(Pos_poly))
    @test_throws MethodError polyhedron(vertices(Pos_poly)) #TODO
    @test_throws ArgumentError convex_hull(rays(Pos_poly))
    @test_throws ArgumentError convex_hull(rays(Pos_poly), [-1 -1 -1])
    @test_throws ArgumentError convex_hull([0 0 0], vertices(Pos_poly))
    @test_throws ArgumentError positive_hull(vertices(Pos_poly))
    @test_throws ArgumentError convex_hull(collect(rays(Pos_poly)))
    @test_throws ArgumentError convex_hull(vertices(Pos_poly), collect(vertices(Pos_poly)))
    @test_throws ArgumentError positive_hull(collect(vertices(Pos_poly)))

    @test_throws ArgumentError incidence_matrix(lineality_space(Pos_poly))
    IM = incidence_matrix([[1]])
    lincone = positive_hull([1 0 0], [0 1 0])

    @test positive_hull(rays_modulo_lineality(lincone)...) == lincone
    @test ambient_dim(polyhedral_fan(IM, rays_modulo_lineality(lincone)...)) == 3
  end
end # of @testset "OscarPolytope"
