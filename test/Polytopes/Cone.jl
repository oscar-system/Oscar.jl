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
        @test size(point_matrix(hilbert_basis(Cone1))) == (2,2)
        @test nrays(Cone1) == 2
        @test point_matrix(rays(Cone1)) == [1 0; 0 1]
        @test point_matrix(facets(Cone1)) == [1 0; 0 1]

        @test !ispointed(Cone2)
        @test !ispointed(Cone3)
        @test !isfulldimensional(Cone4)
        @test isfulldimensional(Cone2)
        @test lineality_space(Cone2)[1] == [0, 1, 0]
        @test Cone2 == Cone3
        @test Cone4 != Cone2
        @test dim(Cone4) == 2
        @test dim(Cone2) == 3
        @test ambient_dim(Cone2) == 3
        @test length(lineality_space(Cone2)) == 1
        @test length(collect(rays(Cone4))) == 2
        @test codim(Cone4) == 1
        @test codim(Cone3) == 0
    end
end
