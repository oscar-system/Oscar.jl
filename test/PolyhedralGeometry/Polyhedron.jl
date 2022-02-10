#TODO: include more examples with nontrivial lineality space

@testset "Polyhedron{$T}" for T in [fmpq, nf_scalar]

    pts = [1 0 0; 0 0 1]'
    @test convex_hull(pts; scalar_type = T) isa Polyhedron{T}
    Q0 = convex_hull(pts; scalar_type = T)
    @test convex_hull(pts; scalar_type = T, non_redundant = true) == Q0
    Q1 = convex_hull(pts, [1 1]; scalar_type = T)
    Q2 = convex_hull(pts, [1 1], [1 1]; scalar_type = T)
    square = cube(T, 2)
    C1 = cube(T, 2, 0, 1)
    Pos = Polyhedron{T}([-1 0 0; 0 -1 0; 0 0 -1], [0,0,0])
    L = Polyhedron{T}([-1 0 0; 0 -1 0], [0,0])
    point = convex_hull([0 1 0]; scalar_type = T)
    s = simplex(T, 2)

    @testset "core functionality" begin
        @test nvertices(Q0) == 3
        @test nvertices.(faces(Q0,1)) == [2,2,2]
        if T == fmpq
            @test lattice_points(Q0) isa SubObjectIterator{PointVector{fmpz}}
            @test point_matrix(lattice_points(Q0)) == matrix(ZZ, [0 0; 0 1; 1 0])
            @test matrix(ZZ, lattice_points(Q0)) == matrix(ZZ, [0 0; 0 1; 1 0])
            @test length(lattice_points(Q0)) == 3
            @test lattice_points(Q0) == [[0, 0], [0, 1], [1, 0]]
            @test interior_lattice_points(square) isa SubObjectIterator{PointVector{fmpz}}
            @test point_matrix(interior_lattice_points(square)) == matrix(ZZ, [0 0])
            @test matrix(ZZ, interior_lattice_points(square)) == matrix(ZZ,[0 0])
            @test length(interior_lattice_points(square)) == 1
            @test interior_lattice_points(square) == [[0, 0]]
            @test boundary_lattice_points(square) isa SubObjectIterator{PointVector{fmpz}}
            @test point_matrix(boundary_lattice_points(square)) == matrix(ZZ, [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1])
            @test length(boundary_lattice_points(square)) == 8
            @test boundary_lattice_points(square) == [[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 1], [1, -1], [1, 0], [1, 1]]
            @test issmooth(Q0)
            @test isnormal(Q0)
        end
        @test isfeasible(Q0)
        @test isbounded(Q0)
        @test isfulldimensional(Q0)
        @test f_vector(Q0) == [3,3]
        @test intersect(Q0, Q0) == Q0
        @test minkowski_sum(Q0, Q0) == convex_hull(2 * pts; scalar_type = T)
        @test Q0+Q0 == minkowski_sum(Q0, Q0)
        @test f_vector(Pos) == [1,3,3]
        @test f_vector(L) == [0, 1, 2]
        @test codim(square) == 0
        @test codim(point) == 3
        @test !isfulldimensional(point)
        @test recession_cone(Pos) isa Cone{T}
        @test nrays(recession_cone(Pos)) == 3
        @test vertices(PointVector{T}, point) isa SubObjectIterator{PointVector{T}}
        @test vertices(PointVector, point) isa SubObjectIterator{PointVector{T}}
        @test vertices(point) isa SubObjectIterator{PointVector{T}}
        if T == fmpq
            @test point_matrix(vertices(2*point)) == matrix(QQ, [0 2 0])
            @test point_matrix(vertices([0,1,0] + point)) == matrix(QQ, [0 2 0])
        else
            @test point_matrix(vertices(2*point)) == [0 2 0]
            @test point_matrix(vertices([0,1,0] + point)) == [0 2 0]
        end
        @test rays(RayVector{T}, Pos) isa SubObjectIterator{RayVector{T}}
        @test rays(RayVector, Pos) isa SubObjectIterator{RayVector{T}}
        @test rays(Pos) isa SubObjectIterator{RayVector{T}}
        @test length(rays(Pos)) == 3
        if T == fmpq
            @test vector_matrix(rays(Pos)) == matrix(QQ, [1 0 0; 0 1 0; 0 0 1])
            @test rays(Pos) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        else
            @test vector_matrix(rays(Pos)) == [0 1 0; 1 0 0; 0 0 1]
            @test rays(Pos) == [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
        end
        @test lineality_space(L) isa SubObjectIterator{RayVector{T}}
        if T == fmpq
            @test generator_matrix(lineality_space(L)) == matrix(QQ, [0 0 1])
            @test matrix(ZZ, lineality_space(L)) == matrix(ZZ, [0 0 1])
        else
            @test generator_matrix(lineality_space(L)) == [0 0 1]
        end
        @test length(lineality_space(L)) == 1
        @test lineality_space(L) == [[0, 0, 1]]
        @test faces(square, 1) isa SubObjectIterator{Polyhedron{T}}
        @test length(faces(square, 1)) == 4
        @test faces(square, 1) == convex_hull.([[-1 -1; -1 1], [1 -1; 1 1], [-1 -1; 1 -1], [-1 1; 1 1]]; scalar_type = T)
        @test vertex_indices(faces(square, 1)) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
        @test ray_indices(faces(square, 1)) == IncidenceMatrix(4, 0)
        @test faces(Pos, 1) isa SubObjectIterator{Polyhedron{T}}
        @test length(faces(Pos, 1)) == 3
        if T == fmpq
            @test faces(Pos, 1) == convex_hull.([[0 0 0]], [[1 0 0], [0 1 0] , [0 0 1]]; scalar_type = T)
        else
            @test faces(Pos, 1) == convex_hull.([[0 0 0]], [[0 1 0], [1 0 0], [0 0 1]]; scalar_type = T)
        end
        @test vertex_indices(faces(Pos, 1)) == IncidenceMatrix([[1], [1], [1]])
        @test ray_indices(faces(Pos, 1)) == IncidenceMatrix([[1], [2], [3]])
        @test isnothing(faces(Q2, 0))
        v = vertices(minkowski_sum(Q0, square))
        @test length(v) == 5
        @test v == [[2, -1], [2, 1], [-1, -1], [-1, 2], [1, 2]]
        if T == fmpq
            @test point_matrix(v) == matrix(QQ, [2 -1; 2 1; -1 -1; -1 2; 1 2])
        else
            @test point_matrix(v) == [2 -1; 2 1; -1 -1; -1 2; 1 2]
        end
        for S in [AffineHalfspace{T}, Pair{Matrix{T}, T}, Polyhedron{T}]
            @test facets(S, Pos) isa SubObjectIterator{S}
            @test length(facets(S, Pos)) == 3
            if T == fmpq
                @test affine_inequality_matrix(facets(S, Pos)) == matrix(QQ, [0 -1 0 0; 0 0 -1 0; 0 0 0 -1])
                @test halfspace_matrix_pair(facets(S, Pos)).A == matrix(QQ, [-1 0 0; 0 -1 0; 0 0 -1]) && halfspace_matrix_pair(facets(S, Pos)).b == [0, 0, 0]
                @test ray_indices(facets(S, Pos)) == IncidenceMatrix([[2, 3], [1, 3], [1, 2]])
            else
                @test affine_inequality_matrix(facets(S, Pos)) == [0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
                @test halfspace_matrix_pair(facets(S, Pos)).A == [-1 0 0; 0 -1 0; 0 0 -1] && halfspace_matrix_pair(facets(S, Pos)).b == [0, 0, 0]
                @test ray_indices(facets(S, Pos)) == IncidenceMatrix([[1, 3], [2, 3], [1, 2]])
            end
            @test vertex_indices(facets(S, Pos)) == IncidenceMatrix([[1], [1], [1]])
            @test facets(S, Pos) == S.([[-1 0 0], [0 -1 0], [0 0 -1]], [0])
        end
        @test facets(Pair, Pos) isa SubObjectIterator{Pair{Matrix{T}, T}}
        @test facets(Pos) isa SubObjectIterator{AffineHalfspace{T}}
        @test facets(Halfspace, Pos) isa SubObjectIterator{AffineHalfspace{T}}
        @test affine_hull(point) isa SubObjectIterator{AffineHyperplane{T}}
        if T == fmpq
            @test affine_equation_matrix(affine_hull(point)) == matrix(QQ, [0 -1 0 0; 1 0 -1 0; 0 0 0 -1])
        else
            @test affine_equation_matrix(affine_hull(point)) == [0 -1 0 0; 1 0 -1 0; 0 0 0 -1]
        end
        @test Oscar.affine_matrix_for_polymake(affine_hull(point)) == [0 -1 0 0; 1 0 -1 0; 0 0 0 -1]
        @test length(affine_hull(point)) == 3
        # TODO: restrict comparison to same scalar?
        @test affine_hull(point) == [Hyperplane([1 0 0], 0), Hyperplane([0 1 0], 1), Hyperplane([0 0 1], 0)]
        @test nfacets(square) == 4
        @test lineality_dim(Q0) == 0
        @test nrays(Q1) == 1
        @test lineality_dim(Q2) == 1
    end

    @testset "linear programs" begin
        LP1 = LinearProgram(square,[1,3])
        LP2 = LinearProgram(square,[2,2]; k=3, convention = :min)
        LP3 = LinearProgram(Pos,[1,2,3])
        @test LP1 isa LinearProgram{T}
        @test LP2 isa LinearProgram{T}
        @test LP3 isa LinearProgram{T}

        @test solve_lp(LP1)==(4,[1,1])
        @test solve_lp(LP2)==(-1,[-1,-1])
        if T == fmpq
            str = ""
        else
            str = "pm::QuadraticExtension<pm::Rational>\n"
        end
        @test string(solve_lp(LP3))==string("(", str, "inf, nothing)")
    end

    @testset "volume" begin
        @test volume(square) isa T
        @test volume(square) == 4
        @test normalized_volume(square) isa T
        @test normalized_volume(square) == 8
        @test normalized_volume(s) == 1
    end

    @testset "standard_constructions" begin
        nc = normal_cone(square, 1)
        @test nc isa Cone{T}
        @test rays(nc) == [[1, 0], [0, 1]]
        if T == fmpq
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
            @test Polyhedron(facets(A)) == A
            b1 = birkhoff(3)
            b2 = birkhoff(3, even = true)
            @test nvertices(pyramid(b1)) + 1 == nvertices(bipyramid(b1))
            @test nvertices(b1) == nvertices(b2) * 2
        end
        P = gelfand_tsetlin([3,2,1]; scalar_type = T)
        p = project_full(P)
        @test p isa Polyhedron{T}
        @test volume(P) == 0
        @test volume(p) == 1
        
    end

end
