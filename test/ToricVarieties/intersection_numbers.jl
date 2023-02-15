using Oscar
using Test

@testset "Topological intersection numbers (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]); set_attributes)

    antv2 = NormalToricVariety([[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]], [[1, 2, 3, 4]]; set_attributes)

    v = NormalToricVariety([[1, 0], [0, 1], [-1, -1]], [[1], [2], [3]]; set_attributes)
    
    dP1 = del_pezzo_surface(1; set_attributes)
    c0 = CohomologyClass(dP1, gens(cohomology_ring(dP1))[1])
    
    dP3 = del_pezzo_surface(3; set_attributes)
    (x1, e1, x2, e3, x3, e2) = gens(cohomology_ring(dP3))
    c1 = CohomologyClass(dP3, x1)
    c2 = CohomologyClass(dP3, e1)
    
    product_space = hirzebruch_surface(5; set_attributes) * projective_space(NormalToricVariety, 2; set_attributes)
    
    @testset "Should fail" begin
        R,_ = PolynomialRing(QQ, 3)
        @test_throws ArgumentError cohomology_ring(antv)
        @test_throws ArgumentError chow_ring(antv2)
        @test_throws ArgumentError c1 - c0
        @test_throws ArgumentError c1 + c0
        @test_throws ArgumentError ideal_of_linear_relations(R, dP3)
    end
    
    @testset "Chow ring and volume form of direct product space" begin
        @test ngens(ideal_of_linear_relations(product_space)) == 4
        @test ngens(chow_ring(product_space).I) == 7
        @test is_trivial(volume_form(product_space)) == false
    end
    
    @testset "Chow ring for non-complete but simplicial varieties" begin
        @test is_simplicial(v) == true
        @test is_complete(v) == false
        @test length(gens(chow_ring(v).I)) == 5
    end
    
    @testset "Properties, attributes and arithmetics of cohomology classes" begin
        @test is_trivial(c1) == false
        @test nrows(exponents(c1)) == 3
        @test length(coefficients(c1)) == 3
        @test fmpq(3) * c1 == fmpz(3) * c1
        @test 2 * c1 != fmpz(3) * c2
        @test (c0 == c1) == false
    end
    
    @testset "Intersection numbers on dP3" begin
        @test integrate(CohomologyClass(dP3, e1*e1)) == -1
        @test integrate(CohomologyClass(dP3, e2*e2)) == -1
        @test integrate(CohomologyClass(dP3, e3*e3)) == -1
        @test integrate(CohomologyClass(dP3, x1*x1)) == -1
        @test integrate(CohomologyClass(dP3, x2*x2)) == -1
        @test integrate(CohomologyClass(dP3, x3*x3)) == -1
        @test integrate(c1) == 0
        @test integrate(c1^2+c1-3//4*c1*c1) == -1//4
        @test length(intersection_form(dP3)) == 21
    end
end
