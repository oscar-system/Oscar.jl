const pm = Polymake


@testset "Cone" begin
    pts = [1 0 0; 0 0 1]'
    Cone1=positive_hull(pts)
    R = [1 0 0; 0 0 1]
    L = [0 1 0]
    Cone2 = Cone(R, L)
    Cone3 = Cone(R, L; non_redundant=true)
    Cone4 = positive_hull(R)
    Cone5 = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
    Cone6 = positive_hull([1//3 1//2; 4//5 2])

    @testset "core functionality" begin
        @test ispointed(Cone1)
        @test isfulldimensional(Cone1)
        @test hilbert_basis(Cone1) isa SubObjectIterator{PointVector{fmpz}}
        @test length(hilbert_basis(Cone1)) == 2
        @test hilbert_basis(Cone1)[1] == PointVector{fmpz}([1, 0])
        @test hilbert_basis(Cone1)[2] == PointVector{fmpz}([0, 1])
        @test generator_matrix(hilbert_basis(Cone1)) == matrix(QQ, [1 0; 0 1])
        @test nrays(Cone1) == 2
        @test rays(RayVector{fmpq}, Cone1) isa SubObjectIterator{RayVector{fmpq}}
        @test rays(Cone1) isa SubObjectIterator{RayVector{fmpq}}
        @test rays(RayVector, Cone1) isa SubObjectIterator{RayVector{fmpq}}
        @test vector_matrix(rays(Cone1)) == matrix(QQ, [1 0; 0 1])
        @test matrix(QQ,rays(Cone1)) == matrix(QQ, [1 0; 0 1])
        @test matrix(ZZ,rays(Cone6)) == matrix(ZZ, [2 3; 2 5])
        @test length(rays(Cone1)) == 2
        @test rays(Cone1)[1] == [1, 0]
        @test rays(Cone1)[2] == [0, 1]
        for T in [AffineHalfspace{fmpq}, LinearHalfspace{fmpq}, Cone{fmpq}, Polyhedron{fmpq}]
            @test facets(T, Cone1) isa SubObjectIterator{T}
            @test length(facets(T, Cone1)) == 2
            @test linear_inequality_matrix(facets(T, Cone1)) == matrix(QQ, [-1 0; 0 -1])
            @test Oscar.linear_matrix_for_polymake(facets(T, Cone1)) == [-1 0; 0 -1]
            @test ray_indices(facets(T, Cone1)) == IncidenceMatrix([[2], [1]])
            if T == Cone{fmpq}
                @test facets(T, Cone1)[1] == cone_from_inequalities([-1 0])
                @test facets(T, Cone1)[2] == cone_from_inequalities([0 -1])
            elseif T == LinearHalfspace{fmpq}
                @test facets(T, Cone1)[1] == T([-1, 0])
                @test facets(T, Cone1)[2] == T([0, -1])
            else
                @test facets(T, Cone1)[1] == T([-1 0], 0)
                @test facets(T, Cone1)[2] == T([0 -1], 0)
            end
        end
        @test facets(Halfspace, Cone1) isa SubObjectIterator{LinearHalfspace{fmpq}}
        @test facets(Cone1) isa SubObjectIterator{LinearHalfspace{fmpq}}
        @test linear_span(Cone4) isa SubObjectIterator{LinearHyperplane{fmpq}}
        @test length(linear_span(Cone4)) == 1
        @test linear_span(Cone4)[] == LinearHyperplane{fmpq}([0 1 0])
        @test linear_equation_matrix(linear_span(Cone4)) == matrix(QQ, [0 1 0])

        @test !ispointed(Cone2)
        @test !ispointed(Cone3)
        @test !isfulldimensional(Cone4)
        @test isfulldimensional(Cone2)
        @test Cone2 == Cone3
        @test Cone4 != Cone2
        @test dim(Cone4) == 2
        @test dim(Cone2) == 3
        @test ambient_dim(Cone2) == 3
        @test lineality_space(Cone2) isa SubObjectIterator{RayVector{fmpq}}
        @test generator_matrix(lineality_space(Cone2)) == matrix(QQ, L)
        @test matrix(QQ, lineality_space(Cone2)) == matrix(QQ, L)
        @test matrix(ZZ, lineality_space(Cone2)) == matrix(ZZ, L)
        @test length(lineality_space(Cone2)) == 1
        @test lineality_space(Cone2)[] == RayVector(L[1, :])
        @test vector_matrix(rays(Cone4)) == matrix(QQ, R)
        @test codim(Cone4) == 1
        @test codim(Cone3) == 0
        @test faces(Cone2, 2) isa SubObjectIterator{Cone{fmpq}}
        @test length(faces(Cone2, 2)) == 2
        @test faces(Cone2, 2)[1] == Cone([0 0 1], [0 1 0])
        @test faces(Cone2, 2)[2] == Cone([1 0 0], [0 1 0])
        @test isnothing(faces(Cone2, 1))
        @test ray_indices(faces(Cone2, 2)) == IncidenceMatrix([[2], [1]])
        @test faces(Cone4, 1) isa SubObjectIterator{Cone{fmpq}}
        @test length(faces(Cone4, 1)) == 2
        @test faces(Cone4, 1)[1] == Cone([0 0 1])
        @test faces(Cone4, 1)[2] == Cone([1 0 0])
        @test ray_indices(faces(Cone4, 1)) == IncidenceMatrix([[2], [1]])

        @test f_vector(Cone5) == [4, 4]
        @test f_vector(Cone2) == [0, 2]
        @test lineality_dim(Cone5) == 0
        @test lineality_dim(Cone2) == 1

        @test nfacets(Cone5) == 4
    end

    @testset "constructors" begin
        @test cone_from_inequalities([-1 0 0; 0 0 -1]) == Cone2
        @test cone_from_inequalities([-1 0 0; 0 0 -1]; non_redundant = true) == Cone2
        @test cone_from_inequalities(facets(Cone4), linear_span(Cone4)) == Cone4
        @test cone_from_inequalities(facets(Cone4), linear_span(Cone4); non_redundant = true) == Cone4
    end
end
