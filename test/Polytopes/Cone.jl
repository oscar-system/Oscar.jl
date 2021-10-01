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
        @test hilbert_basis(Cone1) isa VectorIterator{PointVector{Polymake.Integer}}
        @test (hilbert_basis(Cone1).m == [0 1; 1 0]) || (hilbert_basis(Cone1).m == [1 0; 0 1])
        @test nrays(Cone1) == 2
        @test rays(RayVector{Polymake.Rational}, Cone1) isa VectorIterator{RayVector{Polymake.Rational}}
        @test rays(RayVector{Polymake.Rational}, Cone1).m == [1 0; 0 1]
        @test rays(Cone1) isa VectorIterator{RayVector{Polymake.Rational}}
        @test rays(RayVector, Cone1) isa VectorIterator{RayVector{Polymake.Rational}}
        @test facets(Halfspace, Cone1) isa HalfspaceIterator{Halfspace}
        @test facets(Halfspace, Cone1).A == [-1 0; 0 -1]
        @test facets(Cone, Cone1) isa HalfspaceIterator{Cone}
        @test facets(Cone1) isa HalfspaceIterator{Halfspace}
        @test linear_span(Cone4) isa HalfspaceIterator{Hyperplane}
        @test linear_span(Cone4).A == [0 1 0] && linear_span(Cone4).b == [0]

        @test !ispointed(Cone2)
        @test !ispointed(Cone3)
        @test !isfulldimensional(Cone4)
        @test isfulldimensional(Cone2)
        @test lineality_space(Cone2).m == [0 1 0]
        @test Cone2 == Cone3
        @test Cone4 != Cone2
        @test dim(Cone4) == 2
        @test dim(Cone2) == 3
        @test ambient_dim(Cone2) == 3
        @test lineality_space(Cone2) isa VectorIterator{RayVector{Polymake.Rational}}
        @test lineality_space(Cone2).m == L
        @test rays(Cone4).m == R
        @test codim(Cone4) == 1
        @test codim(Cone3) == 0
        @test faces(Cone, Cone2, 2) isa PolyhedronOrConeIterator{Cone}
        @test size(faces(Cone, Cone2, 2)) == (2,)
        @test faces(Cone, Cone2, 2).lineality == [0 1 0]
        @test faces(Cone2, 2) isa PolyhedronOrConeIterator{Cone}
        @test isnothing(faces(Cone2, 1))

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
