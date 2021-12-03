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

    @testset "core functionality" begin
        @test ispointed(Cone1)
        @test isfulldimensional(Cone1)
        @test hilbert_basis(Cone1) isa SubObjectIterator{PointVector{Polymake.Integer}}
        @test length(hilbert_basis(Cone1)) == 2
        @test hilbert_basis(Cone1)[1] == PointVector{Polymake.Integer}([1, 0])
        @test hilbert_basis(Cone1)[2] == PointVector{Polymake.Integer}([0, 1])
        @test generator_matrix(hilbert_basis(Cone1)) == matrix(QQ, [1 0; 0 1])
        @test nrays(Cone1) == 2
        @test rays(RayVector{Polymake.Rational}, Cone1) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test rays(Cone1) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test rays(RayVector, Cone1) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test vector_matrix(rays(Cone1)) == matrix(QQ, Matrix{fmpq}([1 0; 0 1]))
        @test length(rays(Cone1)) == 2
        @test rays(Cone1)[1] == RayVector([1, 0])
        @test rays(Cone1)[2] == RayVector([0, 1])
        for T in [AffineHalfspace, LinearHalfspace, Cone, Polyhedron]
            @test facets(T, Cone1) isa SubObjectIterator{T}
            @test length(facets(T, Cone1)) == 2
            @test linear_inequality_matrix(facets(T, Cone1)) == matrix(QQ, Matrix{fmpq}([-1 0; 0 -1]))
            if T == Cone
                @test facets(T, Cone1)[1] == cone_from_inequalities([-1 0])
                @test facets(T, Cone1)[2] == cone_from_inequalities([0 -1])
            elseif T == LinearHalfspace
                @test facets(T, Cone1)[1] == T([-1, 0])
                @test facets(T, Cone1)[2] == T([0, -1])
            else
                @test facets(T, Cone1)[1] == T([-1 0], 0)
                @test facets(T, Cone1)[2] == T([0 -1], 0)
            end
        end
        @test facets(Halfspace, Cone1) isa SubObjectIterator{LinearHalfspace}
        @test facets(Cone1) isa SubObjectIterator{LinearHalfspace}
        @test linear_span(Cone4) isa SubObjectIterator{LinearHyperplane}
        @test length(linear_span(Cone4)) == 1
        @test linear_span(Cone4)[] == LinearHyperplane([0 1 0])
        @test linear_equation_matrix(linear_span(Cone4)) == matrix(QQ, Matrix{fmpq}([0 1 0]))

        @test !ispointed(Cone2)
        @test !ispointed(Cone3)
        @test !isfulldimensional(Cone4)
        @test isfulldimensional(Cone2)
        @test Cone2 == Cone3
        @test Cone4 != Cone2
        @test dim(Cone4) == 2
        @test dim(Cone2) == 3
        @test ambient_dim(Cone2) == 3
        @test lineality_space(Cone2) isa SubObjectIterator{RayVector{Polymake.Rational}}
        @test generator_matrix(lineality_space(Cone2)) == matrix(QQ, Matrix{fmpq}(L))
        @test length(lineality_space(Cone2)) == 1
        @test lineality_space(Cone2)[] == RayVector(L[1, :])
        @test vector_matrix(rays(Cone4)) == matrix(QQ, Matrix{fmpq}(R))
        @test codim(Cone4) == 1
        @test codim(Cone3) == 0
        @test faces(Cone2, 2) isa SubObjectIterator{Cone}
        @test length(faces(Cone2, 2)) == 2
        @test faces(Cone2, 2)[1] == Cone([1 0 0], [0 1 0])
        @test faces(Cone2, 2)[2] == Cone([0 0 1], [0 1 0])
        @test isnothing(faces(Cone2, 1))
        @test ray_incidences(faces(Cone2, 2)) == IncidenceMatrix([[1], [2]])

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
