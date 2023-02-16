using Oscar
using Test

@testset "(Weighted) projective spaces (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    P2 = projective_space(NormalToricVariety, 2; set_attributes)
    WPS1 = weighted_projective_space(NormalToricVariety, [2, 3, 1]; set_attributes)
    WPS2 = weighted_projective_space(NormalToricVariety, [1, 2, 3]; set_attributes)
    WPS3 = weighted_projective_space(NormalToricVariety, [3, 1, 2]; set_attributes)
    WPS4 = weighted_projective_space(NormalToricVariety, [3, 45, 68, 7, 15]; set_attributes)
    WPS5 = weighted_projective_space(NormalToricVariety, [3, 1, 2, 6, 5, 1, 1, 8]; set_attributes)
    
    @testset "Basic properties of P2" begin
        @test is_normal(P2) == true
        @test is_affine(P2) == false
        @test is_projective(P2) == true
        @test is_projective_space(P2) == true
        @test is_smooth(P2) == true
        @test is_complete(P2) == true
        @test is_orbifold(P2) == true
        @test is_simplicial(P2) == true
        @test has_torusfactor(P2) == false
    end
    
    @testset "Basic attributes of P2" begin
        @test betti_number(P2, 0) == 1
        @test betti_number(P2, 1) == 0
        @test betti_number(P2, 2) == 1
        @test betti_number(P2, 3) == 0
        @test betti_number(P2, 4) == 1
        @test ngens(cox_ring(P2)) == 3
        @test length(stanley_reisner_ideal(P2).gens) == 1
        @test length(irrelevant_ideal(P2).gens) == 3
        @test transpose(matrix(ZZ,rays(P2))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(P2))
        @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(P2)) == character_lattice(P2)
        @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(P2)) == torusinvariant_weil_divisor_group(P2)
        @test domain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(P2)) == cartier_divisor_group(P2)
        @test codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(P2)) == picard_group(P2)
        @test domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(P2)) == cartier_divisor_group(P2)
        @test codomain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(P2)) == torusinvariant_weil_divisor_group(P2)
        @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(P2)) == matrix(ZZ, [[1], [1], [1]])
        @test coordinate_names(P2) == ["x1", "x2", "x3"]
    end
    
    @testset "Weighted projective space" begin
        @test is_smooth(WPS1) == false
        @test ngens(cox_ring(WPS1)) == 3
        @test dim(WPS1) == 2
        @test iszero(matrix(ZZ, [[2, 3, 1]]) * matrix(ZZ, rays(WPS1)))
        @test is_smooth(WPS2) == false
        @test ngens(cox_ring(WPS2)) == 3
        @test dim(WPS2) == 2
        @test iszero(matrix(ZZ, [[1, 2, 3]]) * matrix(ZZ, rays(WPS2)))
        @test is_smooth(WPS2) == false
        @test ngens(cox_ring(WPS2)) == 3
        @test dim(WPS3) == 2
        @test iszero(matrix(ZZ, [[3, 1, 2]]) * matrix(ZZ, rays(WPS3)))
        @test is_smooth(WPS4) == false
        @test ngens(cox_ring(WPS4)) == 5
        @test dim(WPS4) == 4
        @test iszero(matrix(ZZ, [[3, 45, 68, 7, 15]]) * matrix(ZZ, rays(WPS4)))
        @test is_smooth(WPS5) == false
        @test ngens(cox_ring(WPS5)) == 8
        @test dim(WPS5) == 7
        @test iszero(matrix(ZZ, [[3, 1, 2, 6, 5, 1, 1, 8]]) * matrix(ZZ, rays(WPS5)))
    end
end
