const pm = Polymake

@testset "Cone{$T}" for T in [QQFieldElem, nf_elem]
    
    pts = [1 0; 0 0; 0 1]
    Cone1=positive_hull(T, pts)
    R = [1 0 0; 0 0 1]
    L = [0 1 0]
    Cone2 = positive_hull(T, R, L)
    Cone3 = positive_hull(T, R, L; non_redundant=true)
    Cone4 = positive_hull(T, R)
    Cone5 = positive_hull(T, [1 0 0; 1 1 0; 1 1 1; 1 0 1])
    Cone6 = positive_hull(T, [1//3 1//2; 4//5 2])
    Cone7 = positive_hull(T, [0 1])
    Cone8 = positive_hull(T, [1 1; 1 -1])

    @testset "core functionality" begin
        @test is_pointed(Cone1)
        @test issubset(Cone7, Cone1)
        @test !issubset(Cone1, Cone7)
        @test [1, 0] in Cone1
        @test !([-1, -1] in Cone1)
        if T == QQFieldElem
          @test !is_smooth(Cone2)
          @test is_smooth(Cone7)
          @test !is_smooth(Cone8)
        end
        @test is_simplicial(Cone7)
        @test !is_simplicial(Cone5)
        @test is_fulldimensional(Cone1)
        if T == QQFieldElem
            @test hilbert_basis(Cone1) isa SubObjectIterator{PointVector{ZZRingElem}}
            @test length(hilbert_basis(Cone1)) == 2
            @test hilbert_basis(Cone1) == [[1, 0], [0, 1]]
            @test generator_matrix(hilbert_basis(Cone1)) == matrix(QQ, [1 0; 0 1])
        end
        @test nrays(Cone1) == 2
        @test rays(RayVector{T}, Cone1) isa SubObjectIterator{RayVector{T}}
        @test rays(Cone1) isa SubObjectIterator{RayVector{T}}
        @test rays(RayVector, Cone1) isa SubObjectIterator{RayVector{T}}
        if T == QQFieldElem
            @test vector_matrix(rays(Cone1)) == matrix(QQ, [1 0; 0 1])
            @test matrix(QQ,rays(Cone1)) == matrix(QQ, [1 0; 0 1])
            @test matrix(ZZ,rays(Cone6)) == matrix(ZZ, [2 3; 2 5])
        else
            @test vector_matrix(rays(Cone1)) == [1 0; 0 1]
        end
        @test length(rays(Cone1)) == 2
        @test rays(Cone1) == [[1, 0], [0, 1]]
        for S in [AffineHalfspace{T}, LinearHalfspace{T}, Cone{T}, Polyhedron{T}]
            @test facets(S, Cone1) isa SubObjectIterator{S}
            @test length(facets(S, Cone1)) == 2
            if T == QQFieldElem
                @test linear_inequality_matrix(facets(S, Cone1)) == matrix(QQ, [-1 0; 0 -1])
                @test Oscar.linear_matrix_for_polymake(facets(S, Cone1)) == [-1 0; 0 -1]
                @test ray_indices(facets(S, Cone1)) == IncidenceMatrix([[2], [1]])
                @test IncidenceMatrix(facets(S, Cone1)) == IncidenceMatrix([[2], [1]])
                if S == Cone{T}
                    @test facets(S, Cone1) == cone_from_inequalities.([[-1 0], [0 -1]])
                elseif S == LinearHalfspace{T}
                    @test facets(S, Cone1) == S.([[-1, 0], [0, -1]])
                elseif S == AffineHalfspace{T}
                    @test facets(S, Cone1) == S.([[-1 0], [0 -1]], [0])
                else
                    @test facets(S, Cone1) == polyhedron.(T, [[-1 0], [0 -1]], [0])
                end
            else
                @test linear_inequality_matrix(facets(S, Cone1)) == [0 -1; -1 0]
                @test Oscar.linear_matrix_for_polymake(facets(S, Cone1)) == [0 -1; -1 0]
                @test ray_indices(facets(S, Cone1)) == IncidenceMatrix([[1], [2]])
                @test IncidenceMatrix(facets(S, Cone1)) == IncidenceMatrix([[1], [2]])
                if S == Cone{T}
                    @test facets(S, Cone1) == cone_from_inequalities.(T, [[0 -1], [-1 0]])
                elseif S == LinearHalfspace{T}
                    @test facets(S, Cone1) == S.([[0, -1], [-1, 0]])
                elseif S == AffineHalfspace{T}
                    @test facets(S, Cone1) == S.([[0 -1], [-1 0]], [0])
                else
                    @test facets(S, Cone1) == polyhedron.(T, [[0 -1], [-1 0]], [0])
                end
            end
        end
        @test facets(IncidenceMatrix, Cone1) == IncidenceMatrix(T == QQFieldElem ? [[2], [1]] : [[1], [2]])
        @test facets(Halfspace, Cone1) isa SubObjectIterator{LinearHalfspace{T}}
        @test facets(Cone1) isa SubObjectIterator{LinearHalfspace{T}}
        @test linear_span(Cone4) isa SubObjectIterator{LinearHyperplane{T}}
        @test length(linear_span(Cone4)) == 1
        @test linear_span(Cone4) == [LinearHyperplane{QQFieldElem}([0 1 0])]
        if T == QQFieldElem
            @test linear_equation_matrix(linear_span(Cone4)) == matrix(QQ, [0 1 0])
        else
            @test linear_equation_matrix(linear_span(Cone4)) == [0 1 0]
        end

        @test !is_pointed(Cone2)
        @test !is_pointed(Cone3)
        @test !is_fulldimensional(Cone4)
        @test is_fulldimensional(Cone2)
        @test Cone2 == Cone3
        @test Cone4 != Cone2
        @test dim(Cone4) == 2
        @test dim(Cone2) == 3
        @test ambient_dim(Cone2) == 3
        @test lineality_space(Cone2) isa SubObjectIterator{RayVector{T}}
        if T == QQFieldElem
            @test generator_matrix(lineality_space(Cone2)) == matrix(QQ, L)
            @test matrix(QQ, lineality_space(Cone2)) == matrix(QQ, L)
            @test matrix(ZZ, lineality_space(Cone2)) == matrix(ZZ, L)
        else
            @test generator_matrix(lineality_space(Cone2)) == L
        end
        @test length(lineality_space(Cone2)) == 1
        @test lineality_space(Cone2) == [L[1, :]]
        if T == QQFieldElem
            @test vector_matrix(rays(Cone4)) == matrix(QQ, R)
        else
            @test vector_matrix(rays(Cone4)) == R
        end
        @test codim(Cone4) == 1
        @test codim(Cone3) == 0
        @test faces(Cone2, 2) isa SubObjectIterator{Cone{T}}
        @test length(faces(Cone2, 2)) == 2
        @test faces(Cone4, 1) isa SubObjectIterator{Cone{T}}
        @test length(faces(Cone4, 1)) == 2
        if T == QQFieldElem
            @test faces(Cone2, 2) == positive_hull.(T, [[0 0 1], [1 0 0]], [[0 1 0]])
            @test ray_indices(faces(Cone2, 2)) == IncidenceMatrix([[2], [1]])
            @test IncidenceMatrix(faces(Cone2, 2)) == IncidenceMatrix([[2], [1]])
            @test faces(IncidenceMatrix, Cone2, 2) == IncidenceMatrix([[2], [1]])
            @test faces(Cone4, 1) == positive_hull.(T, [[0 0 1], [1 0 0]])
            @test ray_indices(faces(Cone4, 1)) == IncidenceMatrix([[2], [1]])
            @test IncidenceMatrix(faces(Cone4, 1)) == IncidenceMatrix([[2], [1]])
            @test faces(IncidenceMatrix, Cone4, 1) == IncidenceMatrix([[2], [1]])
        else
            @test faces(Cone2, 2) == positive_hull.(T, [[1 0 0], [0 0 1]], [[0 1 0]])
            @test ray_indices(faces(Cone2, 2)) == IncidenceMatrix([[1], [2]])
            @test IncidenceMatrix(faces(Cone2, 2)) == IncidenceMatrix([[1], [2]])
            @test faces(IncidenceMatrix, Cone2, 2) == IncidenceMatrix([[1], [2]])
            @test faces(Cone4, 1) == positive_hull.(T, [[1 0 0], [0 0 1]])
            @test ray_indices(faces(Cone4, 1)) == IncidenceMatrix([[1], [2]])
            @test IncidenceMatrix(faces(Cone4, 1)) == IncidenceMatrix([[1], [2]])
            @test faces(IncidenceMatrix, Cone4, 1) == IncidenceMatrix([[1], [2]])
        end
        @test IncidenceMatrix(faces(Cone5, 1)) == IncidenceMatrix([[1], [2], [3], [4]])
        @test isnothing(faces(Cone2, 1))

        @test f_vector(Cone5) == [4, 4]
        @test f_vector(Cone2) == [0, 2]
        @test lineality_dim(Cone5) == 0
        @test lineality_dim(Cone2) == 1

        @test nfacets(Cone5) == 4
    end

    @testset "constructors" begin
        @test cone_from_inequalities(T, [-1 0 0; 0 0 -1]) == Cone2
        @test cone_from_inequalities(T, [-1 0 0; 0 0 -1]; non_redundant = true) == Cone2
        @test cone_from_inequalities(T, facets(Cone4), linear_span(Cone4)) == Cone4
        @test cone_from_inequalities(T, facets(Cone4), linear_span(Cone4); non_redundant = true) == Cone4
        @test cone_from_equations(T, [0 1 0]) == cone_from_inequalities(T, Matrix{Int}(undef, 0, 3), linear_span(Cone4))
    end
end
