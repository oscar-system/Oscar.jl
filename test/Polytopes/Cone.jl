const pm = Polymake


@testset "Cone" begin
    pts = [1 0 0; 0 0 1]'
    Cone1=positive_hull(pts)
    R = [1 0 0; 0 0 1]
    L = [0 1 0]
    Cone2 = Cone(R, L)
    Cone3 = Cone(R, L; non_redundant=true)
    Cone4 = positive_hull(R)

    @testset "core functionality" begin
        @test ispointed(Cone1)
        @test isfulldimensional(Cone1)
        @test hilbert_basis(Cone1) isa VectorIterator{PointVector{Polymake.Integer}}
        @test point_matrix(hilbert_basis(Cone1)) == [0 1; 1 0]
        @test nrays(Cone1) == 2
        @test rays(RayVector{Polymake.Rational}, Cone1) isa VectorIterator{RayVector{Polymake.Rational}}
        @test point_matrix(rays(RayVector{Polymake.Rational}, Cone1)) == [1 0; 0 1]
        @test rays(Cone1) isa VectorIterator{RayVector{Polymake.Rational}}
        @test facets(Halfspace, Cone1) isa HalfspaceIterator{Halfspace}
        @test point_matrix(facets(Halfspace, Cone1)) == [1 0; 0 1]
        @test facets(Cone1) isa HalfspaceIterator{Halfspace}

        @test !ispointed(Cone2)
        @test !ispointed(Cone3)
        @test !isfulldimensional(Cone4)
        @test isfulldimensional(Cone2)
        @test point_matrix(lineality_space(Cone2)) == [0 1 0]
        @test Cone2 == Cone3
        @test Cone4 != Cone2
        @test dim(Cone4) == 2
        @test dim(Cone2) == 3
        @test ambient_dim(Cone2) == 3
        @test lineality_space(Cone2) isa VectorIterator{RayVector{Polymake.Rational}}
        @test point_matrix(lineality_space(Cone2)) == L
        @test point_matrix(rays(Cone4)) == R
        @test codim(Cone4) == 1
        @test codim(Cone3) == 0
        @test faces(Cone, Cone2, 2) isa PolyhedronOrConeIterator{Cone}
        @test size(faces(Cone, Cone2, 2)) == (2,)
        @test faces(Cone, Cone2, 2).lineality == [0 1 0]
        @test faces(Cone2, 2) isa PolyhedronOrConeIterator{Cone}
        @test isnothing(faces(Cone2, 1))
    end
end
