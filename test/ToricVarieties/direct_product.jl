@testset "Direct product" begin

    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    P2 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
    ntv6 = F5 * P2

    @testset "Direct product of toric varieties" begin
        @test is_normal(ntv6) == true
        @test is_affine(ntv6) == false
        @test is_projective(ntv6) == true
        @test is_smooth(ntv6) == true
        @test is_complete(ntv6) == true
        @test has_torusfactor(ntv6) == false
        @test is_orbifold(ntv6) == true
        @test is_simplicial(ntv6) == true
        @test betti_number(ntv6, 0) == 1
        @test betti_number(ntv6, 1) == 0
        @test betti_number(ntv6, 2) == 3
        @test betti_number(ntv6, 3) == 0
        @test betti_number(ntv6, 4) == 4
        @test betti_number(ntv6, 5) == 0
        @test betti_number(ntv6, 6) == 3
        @test betti_number(ntv6, 7) == 0
        @test betti_number(ntv6, 8) == 1
        @test is_projective_space(ntv6) == false
    end
end
