#TODO: include more examples with nontrivial lineality space

@testset "Polyhedron{$T}" for T in [fmpq, nf_elem]

    pts = [1 0 0; 0 0 1]'
    @test convex_hull(T, pts) isa Polyhedron{T}
    Q0 = convex_hull(T, pts)
    @test convex_hull(T, pts; non_redundant = true) == Q0
    Q1 = convex_hull(T, pts, [1 1])
    Q2 = convex_hull(T, pts, [1 1], [1 1])
    square = cube(T, 2)
    C1 = cube(T, 2, 0, 1)
    Pos = Polyhedron{T}([-1 0 0; 0 -1 0; 0 0 -1], [0,0,0])
    L = Polyhedron{T}([-1 0 0; 0 -1 0], [0,0])
    point = convex_hull(T, [0 1 0])
    # this is to make sure the order of some matrices below doesn't change
    Polymake.prefer("beneath_beyond") do
        affine_hull(point)
    end
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
            @test is_smooth(Q0)
            @test is_normal(Q0)
        end
        @test is_feasible(Q0)
        @test is_bounded(Q0)
        @test is_fulldimensional(Q0)
        @test f_vector(Q0) == [3,3]
        @test intersect(Q0, Q0) isa Polyhedron{T}
        @test intersect(Q0, Q0) == Q0
        @test minkowski_sum(Q0, Q0) == convex_hull(T, 2 * pts)
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
        @test faces(square, 1) == convex_hull.(T, [[-1 -1; -1 1], [1 -1; 1 1], [-1 -1; 1 -1], [-1 1; 1 1]])
        @test vertex_indices(faces(square, 1)) == IncidenceMatrix([[1, 3], [2, 4], [1, 2], [3, 4]])
        @test ray_indices(faces(square, 1)) == IncidenceMatrix(4, 0)
        @test faces(Pos, 1) isa SubObjectIterator{Polyhedron{T}}
        @test length(faces(Pos, 1)) == 3
        if T == fmpq
            @test faces(Pos, 1) == convex_hull.(T, [[0 0 0]], [[1 0 0], [0 1 0] , [0 0 1]])
        else
            @test faces(Pos, 1) == convex_hull.(T, [[0 0 0]], [[0 1 0], [1 0 0], [0 0 1]])
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
            if S == Pair{Matrix{nf_elem}, nf_elem}
                @test facets(S, Pos) isa SubObjectIterator{Pair{Matrix{Oscar.nf_scalar}, Oscar.nf_scalar}}
                @test facets(S, Pos) == Pair{Matrix{Oscar.nf_scalar}, Oscar.nf_scalar}.([[-1 0 0], [0 -1 0], [0 0 -1]], [0])
            else
                @test facets(S, Pos) isa SubObjectIterator{S}
                @test facets(S, Pos) == S.([[-1 0 0], [0 -1 0], [0 0 -1]], [0])
            end
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
        end
        if T == nf_elem
            @test facets(Pair, Pos) isa SubObjectIterator{Pair{Matrix{Oscar.nf_scalar}, Oscar.nf_scalar}}
        else
            @test facets(Pair, Pos) isa SubObjectIterator{Pair{Matrix{T}, T}}
        end
        @test facets(Pos) isa SubObjectIterator{AffineHalfspace{T}}
        @test facets(Halfspace, Pos) isa SubObjectIterator{AffineHalfspace{T}}
        @test affine_hull(point) isa SubObjectIterator{AffineHyperplane{T}}
        if T == fmpq
            @test affine_equation_matrix(affine_hull(point)) == matrix(QQ, [0 1 0 0; -1 0 1 0; 0 0 0 1])
        else
            @test affine_equation_matrix(affine_hull(point)) == [0 1 0 0; -1 0 1 0; 0 0 0 1]
        end
        @test Oscar.affine_matrix_for_polymake(affine_hull(point)) == [0 1 0 0; -1 0 1 0; 0 0 0 1]
        @test length(affine_hull(point)) == 3
        # TODO: restrict comparison to same scalar?
        @test affine_hull(point) == [Hyperplane([1 0 0], 0), Hyperplane([0 1 0], 1), Hyperplane([0 0 1], 0)]
        @test nfacets(square) == 4
        @test lineality_dim(Q0) == 0
        @test nrays(Q1) == 1
        @test lineality_dim(Q2) == 1
        @test relative_interior_point(Q0) == [1//3, 1//3]
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
        if T == nf_elem
            @test volume(square) isa Oscar.nf_scalar
            @test normalized_volume(square) isa Oscar.nf_scalar
        else
            @test volume(square) isa T
            @test normalized_volume(square) isa T
        end
        @test volume(square) == 4
        @test normalized_volume(square) == 8
        @test normalized_volume(s) == 1
    end

    @testset "standard_constructions" begin
        @test convex_hull(T, pts, nothing, [1 1]) == Q2
        @test Polyhedron{T}(nothing, ([1 0 0; 0 1 0; 0 0 1], [0, 1, 0])) == point
        nc = normal_cone(square, 1)
        @test nc isa Cone{T}
        @test rays(nc) == [[1, 0], [0, 1]]
        let H = LinearHalfspace{T}([1, 1, 0])
            @test Polyhedron(H) isa Polyhedron{T}
            @test Polyhedron(H) == Polyhedron{T}([1 1 0], 0)
        end
        let H = AffineHalfspace{T}([1, 0, 1], 5)
            @test Polyhedron(H) isa Polyhedron{T}
            @test Polyhedron(H) == Polyhedron{T}([1 0 1], 5)
        end
        let H = LinearHyperplane{T}([0, 1, 1])
            @test Polyhedron(H) isa Polyhedron{T}
            @test Polyhedron(H) == Polyhedron{T}((Polymake.Matrix{Polymake.Rational}(undef, 0, 3), Polymake.Rational[]), ([0 1 1], 0))
        end
        let H = AffineHyperplane{T}([1, 1, 1], 7)
            @test Polyhedron(H) isa Polyhedron{T}
            @test Polyhedron(H) == Polyhedron{T}((Polymake.Matrix{Polymake.Rational}(undef, 0, 3), Polymake.Rational[]), ([1 1 1], 7))
        end
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
            
            P = gelfand_tsetlin([3,2,1])
            p = project_full(P)
            @test p isa Polyhedron{T}
            @test volume(P) == 0
            @test volume(p) == 1
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
            
            D = Polyhedron{T}(Polymake.polytope.dodecahedron())
            @test D isa Polyhedron{T}
            @test Polyhedron(Polymake.polytope.dodecahedron()) == D
            
            @test nvertices(D) == 20
            @test vertices(D) == V
            
            let A = [[a//2+1//2 1 0], [0 a//2+1//2 1], [0 a//2+1//2 -1], [-a//2-1//2 -1 0], [a//2-1//2 0 -1], [-a//2-1//2 1 0], [a//2-1//2 0 1], [-a//2+1//2 0 1], [0 -a//2-1//2 1], [-a//2+1//2 0 -1], [a//2+1//2 -1 0], [0 -a//2-1//2 -1]], b = [a//2 + 1, a//2 + 1, a//2 + 1, a//2 + 1, a//4 + 3//4, a//2 + 1, a//4 + 3//4, a//4 + 3//4, a//2 + 1, a//4 + 3//4, a//2 + 1, a//2 + 1]
                
                for S in [AffineHalfspace{T}, Pair{Matrix{T}, T}, Polyhedron{T}]
                    if S == Pair{Matrix{T}, T}
                        @test facets(S, D) isa SubObjectIterator{Pair{Matrix{Oscar.nf_scalar}, Oscar.nf_scalar}}
                        @test facets(S, D) == [Pair{Matrix{Oscar.nf_scalar}, Oscar.nf_scalar}(A[i], b[i]) for i in 1:12]
                    else
                        @test facets(S, D) isa SubObjectIterator{S}
                        @test facets(S, D) == [S(A[i], b[i]) for i in 1:12]
                    end
                    @test length(facets(S, D)) == 12
                    @test affine_inequality_matrix(facets(S, D)) == hcat(-b, vcat(A...))
                    @test halfspace_matrix_pair(facets(S, D)).A == vcat(A...) && halfspace_matrix_pair(facets(S, D)).b == b
                    @test ray_indices(facets(S, D)) == IncidenceMatrix(12, 0)
                    @test vertex_indices(facets(S, D)) == IncidenceMatrix([[1, 3, 5, 9, 10], [1, 2, 3, 4, 6], [1, 2, 5, 7, 8], [11, 12, 16, 18, 20], [5, 8, 10, 15, 17], [2, 4, 7, 11, 12], [3, 6, 9, 13, 14], [4, 6, 11, 13, 16], [13, 14, 16, 19, 20], [7, 8, 12, 15, 18], [9, 10, 14, 17, 19], [15, 17, 18, 19, 20]])
                end
                
            end
            
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
