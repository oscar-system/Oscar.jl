@testset "Algebraic cycles" begin
    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    D2 = DivisorOfCharacter(F5, [1, 2])
    dP1 = NormalToricVariety([[1, 0], [0, 1], [-1, 0], [-1, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    dP3 = NormalToricVariety([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]])
    l4 = canonical_bundle(dP3)
    DC3 = ToricDivisorClass(dP3, [4, 3, 2, 1])
    (u1, u2, u3, u4) = gens(cohomology_ring(dP1))
    (x1, e1, x2, e3, x3, e2) = gens(cohomology_ring(dP3))
    c3 = CohomologyClass(dP1, u1)
    ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
    antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]))
    (xx1, xx2, yy1, yy2) = gens(cox_ring(ntv));
    sv1 = ClosedSubvarietyOfToricVariety(ntv, [xx1])
    sv2 = ClosedSubvarietyOfToricVariety(ntv, [xx1^2+xx1*xx2+xx2^2, yy2])
    P3 = projective_space(NormalToricVariety, 3)
    sv3 = ClosedSubvarietyOfToricVariety(P3, [gens(cox_ring(P3))[1]^2])
    
    ac1 = RationalEquivalenceClass(D2)
    ac2 = RationalEquivalenceClass(DC3)
    ac3 = RationalEquivalenceClass(l4)
    ac4 = RationalEquivalenceClass(c3)
    ac5 = RationalEquivalenceClass(CohomologyClass(dP3, x3))
    ac6 = RationalEquivalenceClass(ToricLineBundle(ntv, [1, 1]))

    @testset "Should fail" begin

        @test_throws ArgumentError RationalEquivalenceClass(antv, [1, 2, 3])
        @test_throws ArgumentError RationalEquivalenceClass(toric_variety(ac1), [1, 2, 3])
        @test_throws ArgumentError ac1 + ac3
        @test_throws ArgumentError ac1 - ac3
        @test_throws ArgumentError ac1 * ac3
        @test_throws ArgumentError sv1 * ac1
        @test_throws ArgumentError ac2 * sv1
        @test_throws ArgumentError sv1 * sv3
    end

    @testset "Arithmetic" begin
        @test (ac1 == ac2) == false
        @test is_trivial(ac2 - ac3) == false
        @test is_trivial(ac1) == true
    end

    @testset "Basic properties" begin
        @test dim(toric_variety(ac1)) == 2
        @test polynomial(ac1) == 0
        @test parent(representative(ac1)) == cox_ring(toric_variety(ac1))
        @test is_trivial(cohomology_class(3*ac4)) == false
        @test is_trivial(RationalEquivalenceClass(sv2)) == false
        @test length(coefficients(ac1)) == 0
        @test length(components(ac2-ac3)) == 3
        @test is_trivial(ac6*sv1) == false
        @test is_trivial(sv1*ac6) == false
        @test is_trivial(sv1*sv1) == true
        @test length(components(sv3*sv3)) == 1
        @test coefficients(sv3*sv3)[1] == 4
    end

    @testset "Intersection" begin
        @test is_trivial(ac2*ac2) == false
        @test is_trivial(ac2*ac2*ac2) == true
    end
end
