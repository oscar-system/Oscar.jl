using Oscar
using Test

@testset "Algebraic cycles (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]))
    
    P3 = projective_space(NormalToricVariety, 3; set_attributes)
    sv0 = ClosedSubvarietyOfToricVariety(P3, [gens(cox_ring(P3))[1]^2])
    
    ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
    (xx1, xx2, yy1, yy2) = gens(cox_ring(ntv));
    sv1 = ClosedSubvarietyOfToricVariety(ntv, [xx1])
    sv2 = ClosedSubvarietyOfToricVariety(ntv, [xx1^2+xx1*xx2+xx2^2, yy2])
    ac0 = RationalEquivalenceClass(ToricLineBundle(ntv, [1, 1]))
    
    F5 = hirzebruch_surface(5; set_attributes)
    ac1 = RationalEquivalenceClass(DivisorOfCharacter(F5, [1, 2]))
    
    dP1 = del_pezzo_surface(1; set_attributes)
    (u1, u2, u3, u4) = gens(cohomology_ring(dP1))
    ac2 = RationalEquivalenceClass(CohomologyClass(dP1, u1))
    
    dP3 = del_pezzo_surface(3; set_attributes)
    (x1, e1, x2, e3, x3, e2) = gens(cohomology_ring(dP3))
    ac3 = RationalEquivalenceClass(canonical_bundle(dP3))
    ac4 = RationalEquivalenceClass(ToricDivisorClass(dP3, [4, 3, 2, 1]))
    
    @testset "Should fail" begin
        @test_throws ArgumentError RationalEquivalenceClass(antv, [1, 2, 3])
        @test_throws ArgumentError RationalEquivalenceClass(toric_variety(ac1), [1, 2, 3])
        @test_throws ArgumentError ac1 + ac3
        @test_throws ArgumentError ac1 - ac3
        @test_throws ArgumentError ac1 * ac3
        @test_throws ArgumentError sv1 * ac1
        @test_throws ArgumentError ac4 * sv1
        @test_throws ArgumentError sv1 * sv0
    end
    
    @testset "Arithmetic" begin
        @test (ac1 == ac4) == false
        @test is_trivial(ac4 - ac3) == false
        @test is_trivial(ac1) == true
    end
    
    @testset "Basic properties" begin
        @test dim(toric_variety(ac1)) == 2
        @test polynomial(ac1) == 0
        @test parent(representative(ac1)) == cox_ring(toric_variety(ac1))
        @test is_trivial(cohomology_class(3*ac2)) == false
        @test is_trivial(RationalEquivalenceClass(sv2)) == false
        @test length(coefficients(ac1)) == 0
        @test length(components(ac4-ac3)) == 4
        @test is_trivial(ac0*sv1) == false
        @test is_trivial(sv1*ac0) == false
        @test is_trivial(sv1*sv1) == true
        @test length(components(sv0*sv0)) == 1
        @test coefficients(sv0*sv0)[1] == 4
    end
    
    @testset "Intersections" begin
        @test is_trivial(ac4*ac4) == false
        @test is_trivial(ac4*ac4*ac4) == true
    end
end