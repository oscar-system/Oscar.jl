#TODO: include more examples with nontrivial lineality space

NF, sr2 = quadratic_field(2)
ENF, sre2 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
Qx, x = QQ["x"]
K, (a1, a2) = embedded_number_field([x^2 - 2, x^3 - 5], [(0, 2), (0, 2)])

for f in (QQ, ENF)

    T = elem_type(f)

    @testset "Polyhedron{$T}" begin

        pts = [1 0; 0 0; 0 1]
        @test convex_hull(f, pts) isa Polyhedron{T}
        Q0 = convex_hull(f, pts)
        @test convex_hull(f, pts; non_redundant = true) == Q0
        Q1 = convex_hull(f, pts, [1 1])
        Q2 = convex_hull(f, pts, [1 1], [1 1])
        square = cube(f, 2)
        CR = cube(f, 2, 0, 3//2)
        Pos = polyhedron(f, [-1 0 0; 0 -1 0; 0 0 -1], [0,0,0])
        L = polyhedron(f, [-1 0 0; 0 -1 0], [0,0])
        point = convex_hull(f, [0 1 0])
        # this is to make sure the order of some matrices below doesn't change
        Polymake.prefer("beneath_beyond") do
            affine_hull(point)
        end
        s = simplex(f, 2)
        R,x = polynomial_ring(QQ, "x")

        @testset "core functionality" begin
            @test issubset(Q0, Q1)
            @test !issubset(Q1, Q0)
            @test [1, 0] in Q0
            @test !([-1, -1] in Q0)
            @test nvertices(Q0) == 3
            @test nvertices.(faces(Q0,1)) == [2,2,2]
            if T == QQFieldElem
                @test lattice_points(Q0) isa SubObjectIterator{PointVector{ZZRingElem}}
                @test point_matrix(lattice_points(Q0)) == matrix(ZZ, [0 0; 0 1; 1 0])
                @test matrix(ZZ, lattice_points(Q0)) == matrix(ZZ, [0 0; 0 1; 1 0])
                @test length(lattice_points(Q0)) == 3
                @test lattice_points(Q0) == [[0, 0], [0, 1], [1, 0]]
                @test interior_lattice_points(square) isa SubObjectIterator{PointVector{ZZRingElem}}
                @test point_matrix(interior_lattice_points(square)) == matrix(ZZ, [0 0])
                @test matrix(ZZ, interior_lattice_points(square)) == matrix(ZZ,[0 0])
                @test length(interior_lattice_points(square)) == 1
                @test interior_lattice_points(square) == [[0, 0]]
                @test boundary_lattice_points(square) isa SubObjectIterator{PointVector{ZZRingElem}}
                @test point_matrix(boundary_lattice_points(square)) == matrix(ZZ, [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1])
                @test length(boundary_lattice_points(square)) == 8
                @test boundary_lattice_points(square) == [[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 1], [1, -1], [1, 0], [1, 1]]
                @test is_smooth(Q0)
                @test is_normal(Q0)
                @test is_lattice_polytope(Q0)
                @test is_very_ample(square)
                @test is_smooth(square)
                @test ehrhart_polynomial(R, square) == 4*x^2 + 4*x + 1
                @test h_star_polynomial(R, CR) == x^4 + 3*x^3 + 10*x^2 + 3*x + 1
                @test is_normal(square)
                @test_throws ArgumentError ehrhart_polynomial(CR)
                @test_throws ArgumentError is_normal(CR)
                @test_throws ArgumentError is_smooth(Q1)
            end
            @test is_feasible(Q0)
            @test is_bounded(Q0)
            @test is_fulldimensional(Q0)
            @test f_vector(Q0) == [3,3]
            @test intersect(Q0, Q0) isa Polyhedron{T}
            @test intersect(Q0, Q0) == Q0
            Ps = [Q0,Q0,Q0]
            @test intersect(Ps) isa Polyhedron{T}
            @test intersect(Ps...) isa Polyhedron{T}
            @test intersect(Ps) == Q0
            @test intersect(Ps...) == Q0
            @test minkowski_sum(Q0, Q0) == convex_hull(f, 2 * pts)
            @test Q0+Q0 == minkowski_sum(Q0, Q0)
            @test f_vector(Pos) == [1,3,3]
            @test f_vector(L) == [0, 1, 2]
            @test codim(square) == 0
            @test codim(point) == 3
            @test !is_fulldimensional(point)
            @test recession_cone(Pos) isa Cone{T}
            @test nrays(recession_cone(Pos)) == 3
            @test vertices(PointVector{T}, point) isa SubObjectIterator{PointVector{T}}
            @test vertices(PointVector, point) isa SubObjectIterator{PointVector{T}}
            @test vertices(point) isa SubObjectIterator{PointVector{T}}
            @test point_matrix(vertices(2*point)) == matrix(f, [0 2 0])
            @test point_matrix(vertices([0,1,0] + point)) == matrix(f, [0 2 0])
            @test rays(RayVector{T}, Pos) isa SubObjectIterator{RayVector{T}}
            @test rays(RayVector, Pos) isa SubObjectIterator{RayVector{T}}
            @test rays(Pos) isa SubObjectIterator{RayVector{T}}
            @test length(rays(Pos)) == 3
            if T == QQFieldElem
                @test vector_matrix(rays(Pos)) == matrix(QQ, [1 0 0; 0 1 0; 0 0 1])
                @test rays(Pos) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            else
                @test vector_matrix(rays(Pos)) == matrix(f, [0 1 0; 1 0 0; 0 0 1])
                @test rays(Pos) == [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
            end
            @test lineality_space(L) isa SubObjectIterator{RayVector{T}}
            @test generator_matrix(lineality_space(L)) == matrix(f, [0 0 1])
            if T == QQFieldElem
                @test matrix(ZZ, lineality_space(L)) == matrix(ZZ, [0 0 1])
            end
            @test length(lineality_space(L)) == 1
            @test lineality_space(L) == [[0, 0, 1]]
            @test faces(square, 1) isa SubObjectIterator{Polyhedron{T}}
            @test length(faces(square, 1)) == 4
            @test faces(square, 1) == convex_hull.([f], [[-1 -1; -1 1], [1 -1; 1 1], [-1 -1; 1 -1], [-1 1; 1 1]])
            @test vertex_indices(faces(square, 1)) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
            @test ray_indices(faces(square, 1)) == IncidenceMatrix(4, 0)
            @test vertex_and_ray_indices(faces(square, 1)) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
            @test IncidenceMatrix(faces(square, 1)) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
            @test faces(IncidenceMatrix, square, 1) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
            @test facet_indices(vertices(square)) == IncidenceMatrix([[1, 3],[2, 3],[1, 4],[2, 4]])
            @test IncidenceMatrix(vertices(square)) == IncidenceMatrix([[1, 3], [2, 3], [1, 4], [2, 4]])
            @test vertices(IncidenceMatrix, square) == IncidenceMatrix([[1, 3], [2, 3], [1, 4], [2, 4]])
            @test faces(Pos, 1) isa SubObjectIterator{Polyhedron{T}}
            @test length(faces(Pos, 1)) == 3
            if T == QQFieldElem
                @test faces(Pos, 1) == convex_hull.(T, [[0 0 0]], [[1 0 0], [0 1 0] , [0 0 1]])
            else
                @test faces(Pos, 1) == convex_hull.([f], [[0 0 0]], [[0 1 0], [1 0 0], [0 0 1]])
            end
            @test vertex_indices(faces(Pos, 1)) == IncidenceMatrix([[1], [1], [1]])
            @test ray_indices(faces(Pos, 1)) == IncidenceMatrix([[1], [2], [3]])
            @test vertex_and_ray_indices(faces(Pos, 1)) == IncidenceMatrix([[1, 4], [2, 4], [3, 4]])
            @test IncidenceMatrix(faces(Pos, 1)) == IncidenceMatrix([[1, 4], [2, 4], [3, 4]])
            @test faces(IncidenceMatrix, Pos, 1) == IncidenceMatrix([[1, 4], [2, 4], [3, 4]])
            @test isnothing(faces(Q2, 0))
            v = vertices(minkowski_sum(Q0, square))
            @test length(v) == 5
            @test v == [f.([2, -1]), f.([2, 1]), f.([-1, -1]), f.([-1, 2]), f.([1, 2])]
            @test point_matrix(v) == matrix(f, [2 -1; 2 1; -1 -1; -1 2; 1 2])
            for S in [AffineHalfspace{T},
                      Pair{Matrix{T}, T},
                      Polyhedron{T}]
                @test facets(S, Pos) isa SubObjectIterator{S}
                if S == Polyhedron{T}
                    @test facets(S, Pos) == polyhedron.([f], [[-1 0 0], [0 -1 0], [0 0 -1]], [0])
                elseif S == Pair{Matrix{T}, T}
                    @test facets(S, Pos) == S.([f.([-1 0 0]), f.([0 -1 0]), f.([0 0 -1])], [f(0)])
                else
                    @test facets(S, Pos) == affine_halfspace.([f], [[-1 0 0], [0 -1 0], [0 0 -1]], [0])
                end
                @test length(facets(S, Pos)) == 3
                @test affine_inequality_matrix(facets(S, Pos)) == matrix(f, [0 -1 0 0; 0 0 -1 0; 0 0 0 -1])
                if T == QQFieldElem
                    @test halfspace_matrix_pair(facets(S, Pos)).A == matrix(QQ, [-1 0 0; 0 -1 0; 0 0 -1]) && halfspace_matrix_pair(facets(S, Pos)).b == [0, 0, 0]
                    @test ray_indices(facets(S, Pos)) == IncidenceMatrix([[2, 3], [1, 3], [1, 2]])
                    @test vertex_and_ray_indices(facets(S, Pos)) == IncidenceMatrix([[2, 3, 4], [1, 3, 4], [1, 2, 4]])
                    @test IncidenceMatrix(facets(S, Pos)) == IncidenceMatrix([[2, 3, 4], [1, 3, 4], [1, 2, 4]])
                else
                    @test halfspace_matrix_pair(facets(S, Pos)).A == [-1 0 0; 0 -1 0; 0 0 -1] && halfspace_matrix_pair(facets(S, Pos)).b == [0, 0, 0]
                    @test ray_indices(facets(S, Pos)) == IncidenceMatrix([[1, 3], [2, 3], [1, 2]])
                    @test vertex_and_ray_indices(facets(S, Pos)) == IncidenceMatrix([[1, 3, 4], [2, 3, 4], [1, 2, 4]])
                    @test IncidenceMatrix(facets(S, Pos)) == IncidenceMatrix([[1, 3, 4], [2, 3, 4], [1, 2, 4]])
                end
                @test vertex_indices(facets(S, Pos)) == IncidenceMatrix([[1], [1], [1]])
            end
            @test facets(IncidenceMatrix, Pos) == IncidenceMatrix(T == QQFieldElem ? [[2, 3, 4], [1, 3, 4], [1, 2, 4]] : [[1, 3, 4], [2, 3, 4], [1, 2, 4]])
            @test facet_indices(vertices(Pos)) == IncidenceMatrix([[1,2,3]])
            @test IncidenceMatrix(vertices(Pos)) == IncidenceMatrix([[1, 2, 3]])
            @test vertices(IncidenceMatrix, Pos) == IncidenceMatrix([[1, 2, 3]])
            @test  facet_indices(rays(Pos)) == ((T==QQFieldElem) ? IncidenceMatrix([[2, 3],[1, 3],[1, 2]]) : IncidenceMatrix([[1, 3],[2, 3],[1, 2]]))
            @test IncidenceMatrix(rays(Pos)) == ((T==QQFieldElem) ? IncidenceMatrix([[2, 3],[1, 3],[1, 2]]) : IncidenceMatrix([[1, 3],[2, 3],[1, 2]]))
            @test rays(IncidenceMatrix, Pos) == ((T==QQFieldElem) ? IncidenceMatrix([[2, 3],[1, 3],[1, 2]]) : IncidenceMatrix([[1, 3],[2, 3],[1, 2]]))
            @test facets(Pair, Pos) isa SubObjectIterator{Pair{Matrix{T}, T}}
            @test facets(Pos) isa SubObjectIterator{AffineHalfspace{T}}
            @test facets(Halfspace, Pos) isa SubObjectIterator{AffineHalfspace{T}}
            @test affine_hull(point) isa SubObjectIterator{AffineHyperplane{T}}
            @test affine_equation_matrix(affine_hull(point)) == matrix(f, [0 1 0 0; -1 0 1 0; 0 0 0 1])
            @test Oscar.affine_matrix_for_polymake(affine_hull(point)) == [0 1 0 0; -1 0 1 0; 0 0 0 1]
            @test length(affine_hull(point)) == 3
            # TODO: restrict comparison to same scalar?
            @test affine_hull(point) == [hyperplane(f, [1 0 0], 0), hyperplane(f, [0 1 0], 1), hyperplane(f, [0 0 1], 0)]
            @test nfacets(square) == 4
            @test lineality_dim(Q0) == 0
            @test nrays(Q1) == 1
            @test lineality_dim(Q2) == 1
            @test relative_interior_point(Q0) == [1//3, 1//3]
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
                @test polyhedron(H) == polyhedron(f, (Polymake.Matrix{Polymake.Rational}(undef, 0, 3), Polymake.Rational[]), ([0 1 1], 0))
            end
            let H = affine_hyperplane(f, [1, 1, 1], 7)
                @test polyhedron(H) isa Polyhedron{T}
                @test polyhedron(H) == polyhedron(f, (Polymake.Matrix{Polymake.Rational}(undef, 0, 3), Polymake.Rational[]), ([1 1 1], 7))
            end
            if T == QQFieldElem
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
                @test polyhedron(facets(A)) == A
                b1 = birkhoff_polytope(3)
                b2 = birkhoff_polytope(3, even = true)
                @test nvertices(pyramid(b1)) + 1 == nvertices(bipyramid(b1))
                @test nvertices(b1) == nvertices(b2) * 2

                P = gelfand_tsetlin_polytope([3,2,1])
                p = project_full(P)
                @test p isa Polyhedron{T}
                @test volume(P) == 0
                @test volume(p) == 1

                rsph = rand_spherical_polytope(3, 15)
                @test rsph isa Polyhedron{T}
                @test is_simplicial(rsph)
                @test nvertices(rsph) == 15

                rsph_r = rand_spherical_polytope(3, 10; distribution=:exact)
                @test map(x->dot(x,x), vertices(rsph_r)) == ones(QQFieldElem,10)
                @test is_simplicial(rsph_r)
                @test nvertices(rsph_r) == 10

                prec = 20
                rsph_prec = rand_spherical_polytope(3, 20; precision=prec)
                @test rsph_prec isa Polyhedron{T}
                @test is_simplicial(rsph_prec)
                @test nvertices(rsph_prec) == 20
                @test all(map(v->abs(dot(v,v)-1), vertices(rsph_prec)) .< QQFieldElem(2)^-(prec-1))

            end

        end

        if T == nf_elem

            @testset "Dodecahedron" begin

                R, a = quadratic_field(5)

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

                D = Polyhedron{T}(Polymake.polytope.dodecahedron(), R)
                @test D isa Polyhedron{T}
                @test polyhedron(Polymake.polytope.dodecahedron()) == D

                @test nvertices(D) == 20
                @test vertices(D) == V

                let A = [[a//2+1//2 1 0], [0 a//2+1//2 1], [0 a//2+1//2 -1], [-a//2-1//2 -1 0], [a//2-1//2 0 -1], [-a//2-1//2 1 0], [a//2-1//2 0 1], [-a//2+1//2 0 1], [0 -a//2-1//2 1], [-a//2+1//2 0 -1], [a//2+1//2 -1 0], [0 -a//2-1//2 -1]], b = [a//2 + 1, a//2 + 1, a//2 + 1, a//2 + 1, a//4 + 3//4, a//2 + 1, a//4 + 3//4, a//4 + 3//4, a//2 + 1, a//4 + 3//4, a//2 + 1, a//2 + 1]

                    for S in [AffineHalfspace{T},
                              Pair{Matrix{T}, T},
                              Polyhedron{T}]
                        
                        @test facets(S, D) isa SubObjectIterator{S}
                        if S == Pair{Matrix{T}, T}
                            @test facets(S, D) == [Pair(A[i], b[i]) for i in 1:12]
                        elseif S == Polyhedron{T}
                            @test facets(S, D) == [polyhedron(R, A[i], b[i]) for i in 1:12]
                        else
                            @test facets(S, D) == [affine_halfspace(R, A[i], b[i]) for i in 1:12]
                        end
                        @test length(facets(S, D)) == 12
                        @test affine_inequality_matrix(facets(S, D)) == matrix(R, hcat(-b, vcat(A...)))
                        @test halfspace_matrix_pair(facets(S, D)).A == vcat(A...) && halfspace_matrix_pair(facets(S, D)).b == b
                        @test ray_indices(facets(S, D)) == IncidenceMatrix(12, 0)
                        @test vertex_indices(facets(S, D)) == IncidenceMatrix([[1, 3, 5, 9, 10], [1, 2, 3, 4, 6], [1, 2, 5, 7, 8], [11, 12, 16, 18, 20], [5, 8, 10, 15, 17], [2, 4, 7, 11, 12], [3, 6, 9, 13, 14], [4, 6, 11, 13, 16], [13, 14, 16, 19, 20], [7, 8, 12, 15, 18], [9, 10, 14, 17, 19], [15, 17, 18, 19, 20]])
                        @test vertex_and_ray_indices(facets(S, D)) == IncidenceMatrix([[1, 3, 5, 9, 10], [1, 2, 3, 4, 6], [1, 2, 5, 7, 8], [11, 12, 16, 18, 20], [5, 8, 10, 15, 17], [2, 4, 7, 11, 12], [3, 6, 9, 13, 14], [4, 6, 11, 13, 16], [13, 14, 16, 19, 20], [7, 8, 12, 15, 18], [9, 10, 14, 17, 19], [15, 17, 18, 19, 20]])
                        @test IncidenceMatrix(facets(S, D)) == IncidenceMatrix([[1, 3, 5, 9, 10], [1, 2, 3, 4, 6], [1, 2, 5, 7, 8], [11, 12, 16, 18, 20], [5, 8, 10, 15, 17], [2, 4, 7, 11, 12], [3, 6, 9, 13, 14], [4, 6, 11, 13, 16], [13, 14, 16, 19, 20], [7, 8, 12, 15, 18], [9, 10, 14, 17, 19], [15, 17, 18, 19, 20]])
                        end

                end

                @test facets(IncidenceMatrix, D) == IncidenceMatrix([[1, 3, 5, 9, 10], [1, 2, 3, 4, 6], [1, 2, 5, 7, 8], [11, 12, 16, 18, 20], [5, 8, 10, 15, 17], [2, 4, 7, 11, 12], [3, 6, 9, 13, 14], [4, 6, 11, 13, 16], [13, 14, 16, 19, 20], [7, 8, 12, 15, 18], [9, 10, 14, 17, 19], [15, 17, 18, 19, 20]])
                @test facet_indices(rays(D))==IncidenceMatrix(0,12)
                @test IncidenceMatrix(rays(D)) == IncidenceMatrix(0, 12)
                @test rays(IncidenceMatrix, D) == IncidenceMatrix(0, 12)
                @test facet_indices(vertices(D))==IncidenceMatrix([[1, 2, 3],[2, 3, 6],[1, 2, 7],[2, 6, 8],[1, 3, 5],[2, 7, 8],[3, 6, 10],[3, 5, 10],[1, 7, 11],[1, 5, 11],[4, 6, 8],[4, 6, 10],[7, 8, 9],[7, 9, 11],[5, 10, 12],[4, 8, 9],[5, 11, 12],[4, 10, 12],[9, 11, 12],[4, 9, 12]])
                @test IncidenceMatrix(vertices(D)) == IncidenceMatrix([[1, 2, 3],[2, 3, 6],[1, 2, 7],[2, 6, 8],[1, 3, 5],[2, 7, 8],[3, 6, 10],[3, 5, 10],[1, 7, 11],[1, 5, 11],[4, 6, 8],[4, 6, 10],[7, 8, 9],[7, 9, 11],[5, 10, 12],[4, 8, 9],[5, 11, 12],[4, 10, 12],[9, 11, 12],[4, 9, 12]])
                @test vertices(IncidenceMatrix, D) == IncidenceMatrix([[1, 2, 3],[2, 3, 6],[1, 2, 7],[2, 6, 8],[1, 3, 5],[2, 7, 8],[3, 6, 10],[3, 5, 10],[1, 7, 11],[1, 5, 11],[4, 6, 8],[4, 6, 10],[7, 8, 9],[7, 9, 11],[5, 10, 12],[4, 8, 9],[5, 11, 12],[4, 10, 12],[9, 11, 12],[4, 9, 12]])
                @test is_feasible(D)
                @test is_bounded(D)
                @test is_fulldimensional(D)
                @test f_vector(D) == [20, 30, 12]
                @test codim(D) == 0
                @test nrays(recession_cone(D)) == 0
                @test nrays(D) == 0
                @test isempty(rays(D))
                @test lineality_dim(D) == 0
                @test isempty(lineality_space(D))
                @test faces(D, 0) == convex_hull.(T, V)
                @test isempty(affine_hull(D))
                @test relative_interior_point(D) == [0, 0, 0]

                @test platonic_solid("dodecahedron") == D

            end

        end

    end
end
