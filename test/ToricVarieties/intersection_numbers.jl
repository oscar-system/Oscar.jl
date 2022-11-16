@testset "Topological intersection numbers" begin
    dP3 = NormalToricVariety([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]])
    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    dP1 = NormalToricVariety([[1, 0], [0, 1], [-1, 0], [-1, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    P2 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
    ntv6 = F5 * P2
    antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]))

    (u1, u2, u3, u4) = gens(cohomology_ring(dP1))
    (x1, e1, x2, e3, x3, e2) = gens(cohomology_ring(dP3))
    c = CohomologyClass(dP3, x1)
    c2 = CohomologyClass(dP3, e1)
    c3 = CohomologyClass(dP1, u1)


    @testset "Should fail" begin
        R,_ = PolynomialRing(QQ, 3)
        @test_throws ArgumentError cohomology_ring(antv)
        @test_throws ArgumentError chow_ring(antv)
        @test_throws ArgumentError is_trivial(c - c3)
        @test_throws ArgumentError is_trivial(c + c3)
        @test_throws ArgumentError ideal_of_linear_relations(R, dP3)
    end

    @testset "Cohomology ring, Chow ring and volume form of direct product space" begin
        @test ngens(ideal_of_linear_relations(ntv6)) == 4
        @test ngens(chow_ring(ntv6).I) == 7
        @test is_trivial(volume_form(ntv6)) == false
    end

    @testset "Properties, attributes and arithmetics of cohomology classes" begin
        @test is_trivial(c) == false
        @test nrows(exponents(c)) == 3
        @test length(coefficients(c)) == 3
        @test fmpq(3) * c == fmpz(3) * c
        @test 2 * c != fmpz(3) * c2
        @test (c == c3) == false
    end

    @testset "Intersection numbers on dP3" begin
        @test integrate(CohomologyClass(dP3, e1*e1)) == -1
        @test integrate(CohomologyClass(dP3, e2*e2)) == -1
        @test integrate(CohomologyClass(dP3, e3*e3)) == -1
        @test integrate(CohomologyClass(dP3, x1*x1)) == -1
        @test integrate(CohomologyClass(dP3, x2*x2)) == -1
        @test integrate(CohomologyClass(dP3, x3*x3)) == -1
        @test integrate(c) == 0
        @test integrate(c^2+c-3//4*c*c) == -1//4
        @test length(intersection_form(dP3)) == 21
        @test integrate(c^2+c-3//4*c*c) == -1//4
    end
end
