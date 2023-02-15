using Oscar
using Test

@testset "Blowup of toric varieties (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    P2 = projective_space(NormalToricVariety, 2; set_attributes)
    BP2 = blowup_on_ith_minimal_torus_orbit(P2, 1, "e"; set_attributes)
    
    @testset "Basic properties of BP2" begin
        @test is_normal(BP2) == true
        @test is_affine(BP2) == false
        @test is_projective(BP2) == true
        @test is_projective_space(BP2) == false
        @test is_smooth(BP2) == true
        @test is_complete(BP2) == true
        @test is_orbifold(BP2) == true
        @test is_simplicial(BP2) == true
        @test has_torusfactor(BP2) == false
    end
    
    @testset "Basic attributes of BP2" begin
        @test betti_number(BP2, 0) == 1
        @test betti_number(BP2, 1) == 0
        @test betti_number(BP2, 2) == 2
        @test betti_number(BP2, 3) == 0
        @test betti_number(BP2, 4) == 1
        @test euler_characteristic(BP2) == 4
        @test rank(picard_group(BP2)) == 2
    end
end