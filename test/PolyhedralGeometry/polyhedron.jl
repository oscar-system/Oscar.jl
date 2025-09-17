#TODO: include more examples with nontrivial lineality space

@testset "Polyhedron{$T}" for (f, T) in _prepare_scalar_types()
  pts = [1 0; 0 0; 0 1]
  @test convex_hull(f, pts) isa Polyhedron{T}
  Q0 = convex_hull(f, pts)
  @test convex_hull(f, pts; non_redundant=true) == Q0
  Q1 = convex_hull(f, pts, [1 1])
  Q2 = convex_hull(f, pts, [1 1], [1 1])
  square = cube(f, 2)
  CR = cube(f, 2, 0, 3//2)
  Pos = polyhedron(f, [-1 0 0; 0 -1 0; 0 0 -1], [0, 0, 0])
  L = polyhedron(f, [-1 0 0; 0 -1 0], [0, 0])
  full = polyhedron(f, zero_matrix(f, 0, 3), [])
  point = convex_hull(f, [0 1 0])
  s = simplex(f, 2)
  R, x = polynomial_ring(QQ, :x)
  v = T[f(1), f(1)]

  # test that remove_zero_rows uses the correct epsilon
  pm_eps = Polymake._get_global_epsilon()
  @test pm_eps < 0.00005

  pf1 = polyhedron([1.0 0.0; 0.0 1.0; 0.0 0.0], [1.0, 1.0, 0.0])
  # three proper ineq, minus the zero row but plus the 1,0,0... row
  @test nrows(Oscar.pm_object(pf1).INEQUALITIES) == 3
  @test n_vertices(pf1) == 1
  pf2 = polyhedron([1.0 0.0; 0.0 1.0; 0.0 0.0], [1.0, 1.0, pm_eps / 2])
  @test nrows(Oscar.pm_object(pf2).INEQUALITIES) == 3
  @test n_vertices(pf2) == 1
  # non-zero row is not removed
  pf3 = polyhedron([1.0 0.0; 0.0 1.0; 0.0 0.0], [1.0, 1.0, pm_eps * 2])
  @test nrows(Oscar.pm_object(pf3).INEQUALITIES) == 4
  @test n_vertices(pf3) == 1

  @testset "core functionality" begin
    @test matrix(f, rays(Q1)) * v == T[f(2)]
    @test issetequal(matrix(f, vertices(Q1)) * v, T[f(1), f(0), f(1)])
    @test issubset(Q0, Q1)
    @test !issubset(Q1, Q0)
    @test [1, 0] in Q0
    @test !([-1, -1] in Q0)
    @test n_vertices(Q0) == 3
    @test n_vertices.(faces(Q0, 1)) == [2, 2, 2]
    @test lattice_points(Q0) isa SubObjectIterator{PointVector{ZZRingElem}}
    @test length(lattice_points(Q0)) == 3
    @test issetequal(lattice_points(Q0), point_vector.(Ref(ZZ), [[0, 0], [0, 1], [1, 0]]))
    @test point_matrix(lattice_points(Q0)) ==
      _oscar_matrix_from_property(ZZ, lattice_points(Q0))
    @test matrix(ZZ, lattice_points(Q0)) ==
      _oscar_matrix_from_property(ZZ, lattice_points(Q0))
    @test interior_lattice_points(square) isa SubObjectIterator{PointVector{ZZRingElem}}
    @test point_matrix(interior_lattice_points(square)) == matrix(ZZ, [0 0])
    @test matrix(ZZ, interior_lattice_points(square)) == matrix(ZZ, [0 0])
    @test length(interior_lattice_points(square)) == 1
    @test interior_lattice_points(square) == [[0, 0]]
    @test boundary_lattice_points(square) isa SubObjectIterator{PointVector{ZZRingElem}}
    @test length(boundary_lattice_points(square)) == 8
    @test issetequal(boundary_lattice_points(square),
      point_vector.(
        Ref(ZZ), [[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 1], [1, -1], [1, 0], [1, 1]]
      ))
    @test point_matrix(boundary_lattice_points(square)) ==
      _oscar_matrix_from_property(ZZ, boundary_lattice_points(square))
    if T == QQFieldElem
      @test is_smooth(Q0)
      @test is_normal(Q0)
      @test is_lattice_polytope(Q0)
      @test is_very_ample(square)
      @test is_smooth(square)
      @test ehrhart_polynomial(R, square) == 4 * x^2 + 4 * x + 1
      @test h_star_polynomial(R, CR) == x^4 + 3 * x^3 + 10 * x^2 + 3 * x + 1
      @test is_normal(square)
      @test_throws ArgumentError ehrhart_polynomial(CR)
      @test_throws ArgumentError is_normal(CR)
      @test_throws ArgumentError is_smooth(Q1)
    end
    @test is_feasible(Q0)
    @test is_bounded(Q0)
    @test is_fulldimensional(Q0)
    @test f_vector(Q0) == [3, 3]
    @test intersect(Q0, Q0) isa Polyhedron{T}
    @test intersect(Q0, Q0) == Q0
    Ps = [Q0, Q0, Q0]
    @test intersect(Ps) isa Polyhedron{T}
    @test intersect(Ps...) isa Polyhedron{T}
    @test intersect(Ps) == Q0
    @test intersect(Ps...) == Q0
    @test minkowski_sum(Q0, Q0) == convex_hull(f, 2 * pts)
    @test Q0 + Q0 == minkowski_sum(Q0, Q0)
    @test f_vector(Pos) == [1, 3, 3]
    @test f_vector(L) == [0, 1, 2]
    @test codim(square) == 0
    @test codim(point) == 3
    @test !is_fulldimensional(point)
    @test recession_cone(Pos) isa Cone{T}
    @test n_rays(recession_cone(Pos)) == 3
    @test vertices(PointVector{T}, point) isa SubObjectIterator{PointVector{T}}
    @test vertices(PointVector, point) isa SubObjectIterator{PointVector{T}}
    @test vertices(point) isa SubObjectIterator{PointVector{T}}
    @test point_matrix(vertices(2 * point)) == matrix(f, [0 2 0])
    @test point_matrix(vertices([0, 1, 0] + point)) == matrix(f, [0 2 0])
    @test rays(RayVector{T}, Pos) isa SubObjectIterator{RayVector{T}}
    @test rays(RayVector, Pos) isa SubObjectIterator{RayVector{T}}
    @test rays(Pos) isa SubObjectIterator{RayVector{T}}
    @test length(rays(Pos)) == 3
    @test issetequal(rays(Pos), ray_vector.(Ref(f), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    @test vector_matrix(rays(Pos)) == _oscar_matrix_from_property(f, rays(Pos))
    @test lineality_space(L) isa SubObjectIterator{RayVector{T}}
    @test generator_matrix(lineality_space(L)) == matrix(f, [0 0 1])
    if T == QQFieldElem
      @test matrix(ZZ, lineality_space(L)) == matrix(ZZ, [0 0 1])
    end
    @test length(lineality_space(L)) == 1
    @test lineality_space(L) == [[0, 0, 1]]
    @test faces(square, 1) isa SubObjectIterator{Polyhedron{T}}
    @test length(faces(square, 1)) == 4
    @test issetequal(faces(square, 1),
      convex_hull.(Ref(f), [[-1 -1; -1 1], [1 -1; 1 1], [-1 -1; 1 -1], [-1 1; 1 1]]))
    @test _check_im_perm_rows(vertex_indices(faces(square, 1)),
      [[1, 3], [2, 4], [1, 2], [3, 4]])
    @test ray_indices(faces(square, 1)) == incidence_matrix(4, 0)
    @test _check_im_perm_rows(vertex_and_ray_indices(faces(square, 1)),
      [[2, 4], [1, 3], [1, 2], [3, 4]])
    @test _check_im_perm_rows(incidence_matrix(faces(square, 1)),
      [[1, 3], [2, 4], [1, 2], [3, 4]])
    @test _check_im_perm_rows(faces(IncidenceMatrix, square, 1),
      [[1, 3], [2, 4], [1, 2], [3, 4]])
    @test _check_im_perm_rows(facet_indices(vertices(square)),
      [[1, 3], [2, 3], [1, 4], [2, 4]])
    @test _check_im_perm_rows(incidence_matrix(vertices(square)),
      [[1, 3], [2, 3], [1, 4], [2, 4]])
    @test _check_im_perm_rows(vertices(IncidenceMatrix, square),
      [[1, 3], [2, 3], [1, 4], [2, 4]])
    @test faces(Pos, 1) isa SubObjectIterator{Polyhedron{T}}
    @test length(faces(Pos, 1)) == 3
    @test issetequal(
      faces(Pos, 1), convex_hull.(Ref(f), [[0 0 0]], [[1 0 0], [0 1 0], [0 0 1]])
    )
    @test _check_im_perm_rows(vertex_indices(faces(Pos, 1)), [[1], [1], [1]])
    @test _check_im_perm_rows(ray_indices(faces(Pos, 1)), [[1], [2], [3]])
    @test _check_im_perm_rows(
      vertex_and_ray_indices(faces(Pos, 1)), [[1, 4], [2, 4], [3, 4]]
    )
    @test _check_im_perm_rows(incidence_matrix(faces(Pos, 1)), [[1, 4], [2, 4], [3, 4]])
    @test _check_im_perm_rows(faces(IncidenceMatrix, Pos, 1), [[1, 4], [2, 4], [3, 4]])
    @test isnothing(faces(Q2, 0))
    v = vertices(minkowski_sum(Q0, square))
    @test length(v) == 5
    @test issetequal(v, point_vector.(Ref(f), [[2, -1], [2, 1], [-1, -1], [-1, 2], [1, 2]]))
    @test point_matrix(v) == _oscar_matrix_from_property(f, v)
    let S = Pair{Matrix{T},T}
      @test issetequal(
        facets(S, Pos), S.([f.([-1 0 0]), f.([0 -1 0]), f.([0 0 -1])], [f(0)])
      )
    end
    let S = AffineHalfspace{T}
      @test issetequal(
        facets(S, Pos), affine_halfspace.(Ref(f), [[-1 0 0], [0 -1 0], [0 0 -1]], [0])
      )
    end
    for S in [AffineHalfspace{T}, Pair{Matrix{T},T}, Polyhedron{T}]
      @test facets(S, Pos) isa SubObjectIterator{S}
      @test length(facets(S, Pos)) == 3
      @test affine_inequality_matrix(facets(S, Pos)) ==
        matrix(f, [0 -1 0 0; 0 0 -1 0; 0 0 0 -1])
      let hmp = halfspace_matrix_pair(facets(S, Pos))
        @test hmp isa NamedTuple{
          (:A, :b),Tuple{elem_type(matrix_space(f, 0, 0)),Vector{elem_type(f)}}
        }
        @test hmp.A == matrix(f, [-1 0 0; 0 -1 0; 0 0 -1])
        @test hmp.b == f.([0, 0, 0])
      end
      @test _check_im_perm_rows(ray_indices(facets(S, Pos)), [[2, 3], [1, 3], [1, 2]])
      @test _check_im_perm_rows(
        vertex_and_ray_indices(facets(S, Pos)), [[2, 3, 4], [1, 3, 4], [1, 2, 4]]
      )
      @test _check_im_perm_rows(
        incidence_matrix(facets(S, Pos)), [[2, 3, 4], [1, 3, 4], [1, 2, 4]]
      )
      @test _check_im_perm_rows(vertex_indices(facets(S, Pos)), [[1], [1], [1]])
    end
    @test _check_im_perm_rows(facets(IncidenceMatrix, Pos),
      [[2, 3, 4], [1, 3, 4], [1, 2, 4]],
    )
    @test _check_im_perm_rows(facet_indices(vertices(Pos)), [[1, 2, 3]])
    @test _check_im_perm_rows(incidence_matrix(vertices(Pos)), [[1, 2, 3]])
    @test _check_im_perm_rows(vertices(IncidenceMatrix, Pos), [[1, 2, 3]])
    @test _check_im_perm_rows(facet_indices(rays(Pos)), [[1, 3], [2, 3], [1, 2]])
    @test _check_im_perm_rows(incidence_matrix(rays(Pos)), [[2, 3], [1, 3], [1, 2]])
    @test _check_im_perm_rows(rays(IncidenceMatrix, Pos), [[1, 3], [2, 3], [1, 2]])
    @test facets(Pair, Pos) isa SubObjectIterator{Pair{Matrix{T},T}}
    @test facets(Pos) isa SubObjectIterator{AffineHalfspace{T}}
    @test facets(Halfspace, Pos) isa SubObjectIterator{AffineHalfspace{T}}
    @test affine_hull(point) isa SubObjectIterator{AffineHyperplane{T}}
    @test length(affine_hull(point)) == 3
    @test issetequal(affine_hull(point),
      [hyperplane(f, [1 0 0], 0), hyperplane(f, [0 1 0], 1), hyperplane(f, [0 0 1], 0)])
    @test affine_equation_matrix(affine_hull(point)) ==
      _oscar_matrix_from_property(f, affine_hull(point))
    @test Oscar.affine_matrix_for_polymake(affine_hull(point)) ==
      _polymake_matrix_from_property(affine_hull(point))
    @test n_facets(square) == 4
    @test lineality_dim(Q0) == 0
    @test n_rays(Q1) == 1
    @test lineality_dim(Q2) == 1
    @test relative_interior_point(Q0) == [1//3, 1//3]

    # issue #4024
    let op = Polyhedron{T}(
        Polymake.polytope.Polytope{Oscar._scalar_type_to_polymake(T)}(;
          REL_INT_POINT=f.([1//3, -1, -2, -4//3]), CONE_AMBIENT_DIM=4
        ),
        f,
      )
      @test relative_interior_point(op) == point_vector(f, [-3, -6, -4])
    end

    @test facet_sizes(Q0)[1] == 2
    @test sum(facet_sizes(Q1)) == 6
    @test facet_sizes(Q2)[1] == 1
    @test vertex_sizes(Q0)[1] == 2
    @test vertex_sizes(Q1)[1] == 2
    @test length(vertex_sizes(Q2)) == 0

    @test length(unique([cube(2), cube(2), simplex(2), simplex(2)])) == 2

    @test dim(full) == ambient_dim(full)
    @test lineality_dim(full) == 3
    @test length(findall(f -> [1, 0] in f, facets(Hyperplane, Q0))) == 2
    @test length(findall(f -> [1, 0] in f, facets(Halfspace, Q0))) == 3
  end

  @testset "volume" begin
    @test volume(square) isa T
    @test normalized_volume(square) isa T
    @test volume(square) == 4
    @test normalized_volume(square) == 8
    @test normalized_volume(s) == 1
  end

  @testset "standard_constructions" begin
    @test convex_hull(f, pts, nothing, [1 1]) == Q2
    @test polyhedron(f, nothing, ([1 0 0; 0 1 0; 0 0 1], [0, 1, 0])) == point
    nc = normal_cone(square, 1)
    @test nc isa Cone{T}
    @test rays(nc) == [[1, 0], [0, 1]]
    let H = linear_halfspace(f, [1, 1, 0])
      @test polyhedron(H) isa Polyhedron{T}
      @test polyhedron(H) == polyhedron(f, [1 1 0], 0)
    end
    let H = affine_halfspace(f, [1, 0, 1], 5)
      @test polyhedron(H) isa Polyhedron{T}
      @test polyhedron(H) == polyhedron(f, [1 0 1], 5)
    end
    let H = linear_hyperplane(f, [0, 1, 1])
      @test polyhedron(H) isa Polyhedron{T}
      @test polyhedron(H) == polyhedron(
        f,
        (Polymake.Matrix{Polymake.Rational}(undef, 0, 3), Polymake.Rational[]),
        ([0 1 1], 0),
      )
    end
    let H = affine_hyperplane(f, [1, 1, 1], 7)
      @test polyhedron(H) isa Polyhedron{T}
      @test polyhedron(H) == polyhedron(
        f,
        (Polymake.Matrix{Polymake.Rational}(undef, 0, 3), Polymake.Rational[]),
        ([1 1 1], 7),
      )
    end
    if T == QQFieldElem
      @test upper_bound_f_vector(4, 8) == [8, 28, 40, 20]
      @test upper_bound_g_vector(4, 8) == [1, 3, 6]
      @test upper_bound_h_vector(4, 8) == [1, 4, 10, 4, 1]
      A = archimedean_solid("cuboctahedron")
      @test count(F -> n_vertices(F) == 3, faces(A, 2)) == 8

      C = catalan_solid("triakis_tetrahedron")
      @test count(F -> n_vertices(F) == 3, faces(C, 2)) == 12

      @test polyhedron(facets(A)) == A
      b1 = birkhoff_polytope(3)
      b2 = birkhoff_polytope(3; even=true)
      @test n_vertices(pyramid(b1)) + 1 == n_vertices(bipyramid(b1))
      @test n_vertices(b1) == n_vertices(b2) * 2

      GT = gelfand_tsetlin_polytope([3, 2, 1])
      @test ambient_dim(GT) == 6
      pGT = project_full(GT)
      @test pGT isa Polyhedron{T}
      @test volume(pGT) == 1

      rsph = rand_spherical_polytope(3, 15)
      @test rsph isa Polyhedron{T}
      @test is_simplicial(rsph)
      @test n_vertices(rsph) == 15

      rsph_r = rand_spherical_polytope(3, 10; distribution=:exact)
      @test map(x -> dot(x, x), vertices(rsph_r)) == ones(QQFieldElem, 10)
      @test is_simplicial(rsph_r)
      @test n_vertices(rsph_r) == 10

      prec = 20
      rsph_prec = rand_spherical_polytope(3, 20; precision=prec)
      @test rsph_prec isa Polyhedron{T}
      @test is_simplicial(rsph_prec)
      @test n_vertices(rsph_prec) == 20
      @test all(
        map(v -> abs(dot(v, v) - 1), vertices(rsph_prec)) .< QQFieldElem(2)^-(prec - 1)
      )

      @test_throws ArgumentError SIM_body_polytope([])
      @test_throws ArgumentError SIM_body_polytope([1, 2, 3])
      let sim = SIM_body_polytope([3, 2, 1])
        @test sim isa Polyhedron{T}
        @test size(Oscar.pm_object(sim).INEQUALITIES, 1) == 16
      end

      let a = associahedron(4)
        @test a isa Polyhedron{T}
        @test dim(a) == 4
        @test n_facets(a) == 14
      end

      let bmg = binary_markov_graph_polytope([0, 1, 0, 0, 1])
        @test bmg isa Polyhedron{T}
        @test ambient_dim(bmg) == 2
        adj = Oscar.pm_object(bmg).SUM_PRODUCT_GRAPH.ADJACENCY
        @test Polymake.nv(adj) == 12
      end

      @test_throws ArgumentError dwarfed_cube(1)
      let dc = dwarfed_cube(3)
        @test dc isa Polyhedron{T}
        @test dim(dc) == 3
        @test n_facets(dc) == 7
      end

      @test_throws ArgumentError dwarfed_product_polygons(3, 2)
      @test_throws ArgumentError dwarfed_product_polygons(4, 2)
      let dpp = dwarfed_product_polygons(4, 5)
        @test dpp isa Polyhedron{T}
        @test is_bounded(dpp)
      end

      @test_throws ArgumentError lecture_hall_simplex(-1)
      let lhs = lecture_hall_simplex(4)
        @test lhs isa Polyhedron{T}
        @test is_bounded(lhs)
        @test n_vertices(lhs) == 5
      end

      let ez = explicit_zonotope([1 2 3; 4 5 6])
        @test ez isa Polyhedron{T}
        @test size(Oscar.pm_object(ez).POINTS, 1) == 4
      end

      @test_throws ArgumentError cyclic_caratheodory_polytope(1, 1)
      @test_throws ArgumentError cyclic_caratheodory_polytope(2, 1)
      let ccp = cyclic_caratheodory_polytope(2, 3)
        @test ccp isa Polyhedron{T}
        @test is_bounded(ccp)
        @test n_vertices(ccp) == 3
      end

      let fkp = fractional_knapsack_polytope([-1, -1, 2])
        @test fkp isa Polyhedron{T}
        @test is_bounded(fkp)
        @test size(Oscar.pm_object(fkp).INEQUALITIES, 1) == 4
      end

      @test_throws ArgumentError hypersimplex(5, 3)
      let hs = hypersimplex(3, 5)
        @test hs isa Polyhedron{T}
        @test is_bounded(hs)
        @test n_vertices(hs) == 10
        @test length(facets(hs)) == 10
      end

      @test zonotope([0 0 1; 2 2 2; 1 0 2]) isa Polyhedron{T}

      @test_throws ArgumentError goldfarb_cube(0, 0, 0)
      @test_throws ArgumentError goldfarb_cube(3, 1, 0)
      @test_throws ArgumentError goldfarb_cube(3, 1//4, 1//8)
      let goldfarb = goldfarb_cube(3, 1//4, 0) # fuer zuweisungen 
        @test goldfarb isa Polyhedron{T}
        @test ambient_dim(goldfarb) == 3
        @test size(Oscar.pm_object(goldfarb).INEQUALITIES, 1) == 7 #nur falls oscar danach nicht fragen kann
      end

      @test_throws ArgumentError goldfarb_sit_cube(0, 0, 0)
      @test_throws ArgumentError goldfarb_sit_cube(3, 1, 0)
      @test_throws ArgumentError goldfarb_sit_cube(3, 1//4, 5//8)
      let gsc = goldfarb_sit_cube(3, 1//4, 0)
        @test gsc isa Polyhedron{T}
        @test ambient_dim(gsc) == 3
        @test size(Oscar.pm_object(gsc).INEQUALITIES, 1) == 7
      end

      @test_throws ArgumentError hypertruncated_cube(3, 5, 0)
      let htc = hypertruncated_cube(3, 3//2, 3//4)
        @test htc isa Polyhedron{T}
        @test is_bounded(htc)
        @test size(Oscar.pm_object(htc).INEQUALITIES, 1) == 13
      end

      let kcp = k_cyclic_polytope(8, [1, 2])
        @test kcp isa Polyhedron{T}
        @test n_vertices(kcp) == 8
      end

      @test_throws ArgumentError klee_minty_cube(0, 0)
      @test_throws ArgumentError klee_minty_cube(3, 1)
      let kmc = klee_minty_cube(4, 1//32)
        @test kmc isa Polyhedron{T}
        @test is_bounded(kmc)
        @test size(Oscar.pm_object(kmc).INEQUALITIES, 1) == 9
      end

      @test_throws ArgumentError max_GC_rank_polytope(1)
      @test_throws ArgumentError max_GC_rank_polytope(100000)
      let gc = max_GC_rank_polytope(4)
        @test gc isa Polyhedron{T}
        @test ambient_dim(gc) == 4
        @test is_bounded(gc)
      end

      @test_throws ArgumentError n_gon(2)
      let gon = n_gon(4; r=2)
        @test gon isa Polyhedron{QQBarFieldElem}
        @test length(vertices(gon)[1]) == 2
        @test n_vertices(gon) == 4
      end
      let gon = n_gon(20)
        K = coefficient_field(gon)
        @test volume(gon) == 5//2 * (sqrt(K(5)) - 1)
      end

      @test_throws ArgumentError permutahedron(-1)
      let perm = permutahedron(3)
        @test perm isa Polyhedron{T}
        @test length(vertices(perm)[1]) == 4
        @test dim(perm) == 3
      end

      let pile = pile_polytope([2, 2])
        @test pile isa Polyhedron{T}
        @test ambient_dim(pile) == 3
        @test n_vertices(pile) == 9
      end

      @test_throws ArgumentError pitman_stanley_polytope(Vector{Rational}([]))
      let psp = pitman_stanley_polytope([2, 4, 2])
        @test psp isa Polyhedron{T}
        @test ambient_dim(psp) == 3
        @test size(Oscar.pm_object(psp).INEQUALITIES, 1) == 6
      end

      @test_throws ArgumentError pseudo_del_pezzo_polytope(0)
      let pdp = pseudo_del_pezzo_polytope(2)
        @test pdp isa Polyhedron{T}
        @test is_bounded(pdp)
        @test dim(pdp) == 2
      end

      @test_throws ArgumentError rand01_polytope(1, 1)
      @test rand01_polytope(2, 4) isa Polyhedron{T}
      let r_01_p = rand01_polytope(2, 4; seed=47)
        @test r_01_p isa Polyhedron{T}
        @test n_vertices(r_01_p) == 4
      end

      @test_throws ArgumentError rand_box_polytope(1, 0, 1)
      let rbox = rand_box_polytope(3, 8, 1)
        @test rbox isa Polyhedron{T}
        @test ambient_dim(rbox) == 3
        @test size(Oscar.pm_object(rbox).POINTS, 1) == 8
      end
      let rbox = rand_box_polytope(3, 4, 1; seed=456)
        @test rbox isa Polyhedron{T}
        @test sum(Oscar.pm_object(rbox).POINTS) == 6
      end

      let rmetric = rand_metric(3; seed=213)
        @test rmetric isa QQMatrix
      end

      let rmetricint = rand_metric_int(3, 2; seed=213)
        @test rmetricint isa ZZMatrix
      end

      let rnorm = rand_normal_polytope(3, 4; seed=213)
        @test rnorm isa Polyhedron{T}
        @test is_bounded(rnorm)
        @test size(Oscar.pm_object(rnorm).POINTS, 1) == 4
      end

      @test_throws ArgumentError rand_cyclic_polytope(2, 3)
      @test rand_cyclic_polytope(2, 4) isa Polyhedron{T}
      let rcyc = p = rand_cyclic_polytope(3, 8; seed=4)
        @test rcyc isa Polyhedron{T}
        @test n_vertices(rcyc) == 8
      end

      #3 more rand testset

      @test_throws ArgumentError rss_associahedron(1)
      let rss = rss_associahedron(3)
        @test rss isa Polyhedron{T}
        @test typeof(facets(rss)[1]) == AffineHalfspace{T}
        @test ambient_dim(rss) == 3
      end

      @test_throws ArgumentError signed_permutahedron(0)
      @test_throws ArgumentError signed_permutahedron(100000)
      let sph = signed_permutahedron(3)
        @test sph isa Polyhedron{T}
        @test n_vertices(sph) == 48
        @test is_bounded(sph)
      end

      let G = complete_graph(3)
        ssp = stable_set_polytope(G)
        @test ssp isa Polyhedron{T}
        @test is_bounded(ssp)
        @test size(Oscar.pm_object(ssp).INEQUALITIES, 1) == 7
      end

      @test_throws ArgumentError transportation_polytope([2, 3, 1], [4, 5, 3])
      let tp = transportation_polytope([2, 3, 1], [4, 0, 2])
        @test tp isa Polyhedron{T}
        @test ambient_dim(tp) == 9
        @test size(Oscar.pm_object(tp).EQUATIONS, 1) == 6
      end

      @test sum(zonotope_vertices_fukuda_matrix([4 4 2; 2 4 1])) == 0

      @test_throws ArgumentError vertex_figure(cube(3), 10)
      let vf = vertex_figure(platonic_solid("octahedron"), 1)
        @test vf isa Polyhedron{T}
        @test ambient_dim(vf) == 3
        @test sum(facet_sizes(vf)) == 8
      end
    end
  end

  if T == EmbeddedNumFieldElem{AbsSimpleNumFieldElem}
    @testset "Dodecahedron" begin
      D = polyhedron(Polymake.polytope.dodecahedron())
      R = coefficient_field(D)
      NF = number_field(R)
      let isq = Hecke.is_quadratic_type(NF)
        @test isq[1]
        @test isq[2] == 5
      end
      a = R(gens(NF)[])

      V = [[1//2, a//4 + 3//4, 0],
        [-1//2, a//4 + 3//4, 0],
        [a//4 + 1//4, a//4 + 1//4, a//4 + 1//4],
        [-a//4 - 1//4, a//4 + 1//4, a//4 + 1//4],
        [a//4 + 1//4, a//4 + 1//4, -a//4 - 1//4],
        [0, 1//2, a//4 + 3//4],
        [-a//4 - 1//4, a//4 + 1//4, -a//4 - 1//4],
        [0, 1//2, -a//4 - 3//4],
        [a//4 + 3//4, 0, 1//2],
        [a//4 + 3//4, 0, -1//2],
        [-a//4 - 3//4, 0, 1//2],
        [-a//4 - 3//4, 0, -1//2],
        [0, -1//2, a//4 + 3//4],
        [a//4 + 1//4, -a//4 - 1//4, a//4 + 1//4],
        [0, -1//2, -a//4 - 3//4],
        [-a//4 - 1//4, -a//4 - 1//4, a//4 + 1//4],
        [a//4 + 1//4, -a//4 - 1//4, -a//4 - 1//4],
        [-a//4 - 1//4, -a//4 - 1//4, -a//4 - 1//4],
        [1//2, -a//4 - 3//4, 0],
        [-1//2, -a//4 - 3//4, 0]]

      @test D isa Polyhedron{T}

      @test n_vertices(D) == 20
      @test vertices(D) == V

      let A = [
          [a//2 + 1//2 1 0],
          [0 a//2 + 1//2 1],
          [0 a//2 + 1//2 -1],
          [-a//2 - 1//2 -1 0],
          [a//2 - 1//2 0 -1],
          [-a//2 - 1//2 1 0],
          [a//2 - 1//2 0 1],
          [-a//2 + 1//2 0 1],
          [0 -a//2 - 1//2 1],
          [-a//2 + 1//2 0 -1],
          [a//2 + 1//2 -1 0],
          [0 -a//2 - 1//2 -1],
        ],
        b = [
          a//2 + 1,
          a//2 + 1,
          a//2 + 1,
          a//2 + 1,
          a//4 + 3//4,
          a//2 + 1,
          a//4 + 3//4,
          a//4 + 3//4,
          a//2 + 1,
          a//4 + 3//4,
          a//2 + 1,
          a//2 + 1,
        ]

        for S in [AffineHalfspace{T},
          Pair{Matrix{T},T},
          Polyhedron{T}]
          @test facets(S, D) isa SubObjectIterator{S}
          if S == Pair{Matrix{T},T}
            @test facets(S, D) == [Pair(A[i], b[i]) for i in 1:12]
          elseif S == Polyhedron{T}
            @test n_vertices.(facets(S, D)) == repeat([5], 12)
          else
            @test facets(S, D) == [affine_halfspace(R, A[i], b[i]) for i in 1:12]
          end
          @test length(facets(S, D)) == 12
          @test affine_inequality_matrix(facets(S, D)) == matrix(R, hcat(-b, vcat(A...)))
          @test halfspace_matrix_pair(facets(S, D)) == (A=matrix(R, vcat(A...)), b=b)

          @test ray_indices(facets(S, D)) == incidence_matrix(12, 0)
          @test _check_im_perm_rows(
            vertex_indices(facets(S, D)),
            [
              [1, 3, 5, 9, 10],
              [1, 2, 3, 4, 6],
              [1, 2, 5, 7, 8],
              [11, 12, 16, 18, 20],
              [5, 8, 10, 15, 17],
              [2, 4, 7, 11, 12],
              [3, 6, 9, 13, 14],
              [4, 6, 11, 13, 16],
              [13, 14, 16, 19, 20],
              [7, 8, 12, 15, 18],
              [9, 10, 14, 17, 19],
              [15, 17, 18, 19, 20],
            ],
          )
          @test _check_im_perm_rows(
            vertex_and_ray_indices(facets(S, D)),
            [
              [1, 3, 5, 9, 10],
              [1, 2, 3, 4, 6],
              [1, 2, 5, 7, 8],
              [11, 12, 16, 18, 20],
              [5, 8, 10, 15, 17],
              [2, 4, 7, 11, 12],
              [3, 6, 9, 13, 14],
              [4, 6, 11, 13, 16],
              [13, 14, 16, 19, 20],
              [7, 8, 12, 15, 18],
              [9, 10, 14, 17, 19],
              [15, 17, 18, 19, 20],
            ],
          )
          @test _check_im_perm_rows(
            incidence_matrix(facets(S, D)),
            [
              [1, 3, 5, 9, 10],
              [1, 2, 3, 4, 6],
              [1, 2, 5, 7, 8],
              [11, 12, 16, 18, 20],
              [5, 8, 10, 15, 17],
              [2, 4, 7, 11, 12],
              [3, 6, 9, 13, 14],
              [4, 6, 11, 13, 16],
              [13, 14, 16, 19, 20],
              [7, 8, 12, 15, 18],
              [9, 10, 14, 17, 19],
              [15, 17, 18, 19, 20],
            ],
          )
        end
      end

      @test _check_im_perm_rows(
        facets(IncidenceMatrix, D),
        [
          [1, 3, 5, 9, 10],
          [1, 2, 3, 4, 6],
          [1, 2, 5, 7, 8],
          [11, 12, 16, 18, 20],
          [5, 8, 10, 15, 17],
          [2, 4, 7, 11, 12],
          [3, 6, 9, 13, 14],
          [4, 6, 11, 13, 16],
          [13, 14, 16, 19, 20],
          [7, 8, 12, 15, 18],
          [9, 10, 14, 17, 19],
          [15, 17, 18, 19, 20],
        ],
      )
      @test facet_indices(rays(D)) == incidence_matrix(0, 12)
      @test incidence_matrix(rays(D)) == incidence_matrix(0, 12)
      @test rays(IncidenceMatrix, D) == incidence_matrix(0, 12)
      @test _check_im_perm_rows(
        facet_indices(vertices(D)),
        [
          [1, 2, 3],
          [2, 3, 6],
          [1, 2, 7],
          [2, 6, 8],
          [1, 3, 5],
          [2, 7, 8],
          [3, 6, 10],
          [3, 5, 10],
          [1, 7, 11],
          [1, 5, 11],
          [4, 6, 8],
          [4, 6, 10],
          [7, 8, 9],
          [7, 9, 11],
          [5, 10, 12],
          [4, 8, 9],
          [5, 11, 12],
          [4, 10, 12],
          [9, 11, 12],
          [4, 9, 12],
        ],
      )
      @test _check_im_perm_rows(
        incidence_matrix(vertices(D)),
        [
          [1, 2, 3],
          [2, 3, 6],
          [1, 2, 7],
          [2, 6, 8],
          [1, 3, 5],
          [2, 7, 8],
          [3, 6, 10],
          [3, 5, 10],
          [1, 7, 11],
          [1, 5, 11],
          [4, 6, 8],
          [4, 6, 10],
          [7, 8, 9],
          [7, 9, 11],
          [5, 10, 12],
          [4, 8, 9],
          [5, 11, 12],
          [4, 10, 12],
          [9, 11, 12],
          [4, 9, 12],
        ],
      )
      @test _check_im_perm_rows(
        vertices(IncidenceMatrix, D),
        [
          [1, 2, 3],
          [2, 3, 6],
          [1, 2, 7],
          [2, 6, 8],
          [1, 3, 5],
          [2, 7, 8],
          [3, 6, 10],
          [3, 5, 10],
          [1, 7, 11],
          [1, 5, 11],
          [4, 6, 8],
          [4, 6, 10],
          [7, 8, 9],
          [7, 9, 11],
          [5, 10, 12],
          [4, 8, 9],
          [5, 11, 12],
          [4, 10, 12],
          [9, 11, 12],
          [4, 9, 12],
        ],
      )
      @test is_feasible(D)
      @test is_bounded(D)
      @test is_fulldimensional(D)
      @test f_vector(D) == [20, 30, 12]
      @test codim(D) == 0
      @test n_rays(recession_cone(D)) == 0
      @test n_rays(D) == 0
      @test isempty(rays(D))
      @test lineality_dim(D) == 0
      @test isempty(lineality_space(D))
      @test faces(D, 0) == convex_hull.(T, V)
      @test isempty(affine_hull(D))
      @test relative_interior_point(D) == [0, 0, 0]
    end
  end

  @testset "Cone conversion" begin
    C1 = positive_hull(f, [1 0])
    C2 = positive_hull(f, [1 0], [0 1])
    C3 = positive_hull(f, [1 0 0; 1 1 0; 1 1 1; 1 0 1])
    for C in [C1, C2, C3]
      PC = polyhedron(C)
      @test dim(C) == dim(PC)
      @test ambient_dim(C) == ambient_dim(PC)
      @test lineality_dim(C) == lineality_dim(PC)
      @test matrix(f, rays(C)) == matrix(f, rays(PC))
      vdesc = convex_hull(
        minimal_faces(PointVector, PC),
        rays_modulo_lineality(RayVector, PC),
        lineality_space(PC),
      )
      hdesc = polyhedron(facets(PC), affine_hull(PC))
      @test vdesc == hdesc
      @test recession_cone(PC) == C
    end
    pts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
    inc = incidence_matrix([
      [4, 5, 6], [1, 4, 2], [2, 4, 5], [2, 3, 5], [3, 5, 6], [1, 3, 6], [1, 4, 6]
    ])
    SOP = subdivision_of_points(pts, inc)
    SC = secondary_cone(SOP)
    @test polyhedron(SC) isa Polyhedron
  end
end

@testset "Regular solids" begin
  for i in Oscar._johnson_indexes_from_oscar
    j = johnson_solid(i)
    @test j isa Polyhedron{<:EmbeddedNumFieldElem}
    @test Polymake.polytope.isomorphic(
      Oscar.pm_object(j), Polymake.polytope.johnson_solid(i)
    )
  end

  let p = platonic_solid("dodecahedron")
    @test is_platonic_solid(p)
    @test !is_archimedean_solid(p)
    @test !is_johnson_solid(p)
    @test is_vertex_transitive(p)
  end

  let a = archimedean_solid("rhombicuboctahedron")
    @test !is_platonic_solid(a)
    @test is_archimedean_solid(a)
    @test !is_johnson_solid(a)
    @test is_vertex_transitive(a)
    @test !Oscar._is_prismic_or_antiprismic(a)
  end

  for f in (snub_cube, snub_dodecahedron)
    ae = f()
    @test ae isa Polyhedron{<:EmbeddedNumFieldElem}
    @test is_archimedean_solid(ae)
    @test Polymake.polytope.isomorphic(
      Oscar.pm_object(ae), Polymake.polytope.archimedean_solid(string(f))
    )
  end

  for f in (pentagonal_icositetrahedron, pentagonal_hexecontahedron)
    ce = f()
    @test ce isa Polyhedron{<:EmbeddedNumFieldElem}
    @test Polymake.polytope.isomorphic(
      Oscar.pm_object(ce), Polymake.polytope.catalan_solid(string(f))
    )
  end

  let j = johnson_solid(69)
    @test !is_platonic_solid(j)
    @test !is_archimedean_solid(j)
    @test is_johnson_solid(j)
    @test !is_vertex_transitive(j)
  end

  K = algebraic_closure(QQ)
  e = one(K)
  s, c = sinpi(2 * e / 7), cospi(2 * e / 7)
  mat_rot = matrix([c -s; s c])
  mat_sigma1 = matrix(K, [-1 0; 0 1])
  G_mat = matrix_group(mat_rot, mat_sigma1)
  p = K.([0, 1])  # coordinates of vertex 1, expressed over K
  orb = orbit(G_mat, *, p)
  pts = collect(orb)
  len = sqrt(dot(pts[1] - pts[2], pts[1] - pts[2])) / 2
  ngon = convex_hull(pts)
  prism = polyhedron(Polymake.polytope.prism(Oscar.pm_object(ngon), len))
  @test Oscar._is_prismic_or_antiprismic(prism)
  @test !is_archimedean_solid(prism)
  @test !is_johnson_solid(prism)
end
