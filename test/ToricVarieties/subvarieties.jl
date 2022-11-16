@testset "Closed subvarieties" begin
    ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
    P3 = projective_space(NormalToricVariety, 3)
    (x1, x2, y1, y2) = gens(cox_ring(ntv));


    @testset "Should fail" begin
        antv6 = NormalToricVariety([[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]], [[1, 2, 3, 4]])
        @test_throws ArgumentError ClosedSubvarietyOfToricVariety(ntv, [x1 - y1])
        @test_throws ArgumentError ClosedSubvarietyOfToricVariety(antv6, [gens(cox_ring(antv6))[1]])
    end

    @testset "Basic properties" begin
        sv1 = ClosedSubvarietyOfToricVariety(ntv, [x1])
        sv2 = ClosedSubvarietyOfToricVariety(ntv, [x1^2+x1*x2+x2^2, y2])
        sv3 = ClosedSubvarietyOfToricVariety(P3, [gens(cox_ring(P3))[1]^2])
        @test sv1 isa ClosedSubvarietyOfToricVariety
        @test sv2 isa ClosedSubvarietyOfToricVariety
        @test sv3 isa ClosedSubvarietyOfToricVariety
        @test is_empty(sv1) == false
        @test radical(sv1) == defining_ideal(sv1)
        @test dim(toric_variety(sv1)) == 2
    end
end
