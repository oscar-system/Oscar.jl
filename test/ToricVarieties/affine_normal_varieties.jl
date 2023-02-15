using Oscar
using Test

@testset "Affine normal toric varieties (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    antv = AffineNormalToricVariety(Oscar.positive_hull([1 1; -1 1]); set_attributes)
    antv2 = NormalToricVariety(Oscar.positive_hull([1 1; -1 1]); set_attributes)
    antv3 = AffineNormalToricVariety(antv2; set_attributes)
    antv4 = AffineNormalToricVariety(Oscar.positive_hull([1 0]); set_attributes)
    antv5 = affine_space(NormalToricVariety, 2; set_attributes)
    antv6 = NormalToricVariety([[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]], [[1, 2, 3, 4]]; set_attributes)
    
    @testset "Basic properties" begin
        @test is_smooth(antv) == false
        @test is_orbifold(antv) == true
        @test is_projective_space(antv) == false
    end
    
    @testset "Basic attributes" begin
        @test dim(fan(antv)) == 2
        @test dim(cone(antv)) == 2
        @test length(affine_open_covering(antv)) == 1
        @test length(gens(toric_ideal(antv))) == 1
        @test rank(torusinvariant_weil_divisor_group(antv)) == 2
        @test rank(character_lattice(antv)) == 2
        @test rank(domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(antv))) == 2
        @test rank(codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(antv))) == 2
        @test elementary_divisors(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(antv))) == [ 2 ]
        @test elementary_divisors(class_group(antv)) == [ 2 ]
        @test ngens(cox_ring(antv)) == 2
        @test length(torusinvariant_prime_divisors(antv)) == 2
    end
    
    @testset "Cast to affine" begin
        @test is_affine(antv2) == true
        @test ngens(toric_ideal(antv2)) == 1
        @test length(affine_open_covering(antv3)) == 1
        @test dim(cone(antv3)) == 2
    end
    
    @testset "Finalizing a variety with torus-factor" begin
        @test dim_of_torusfactor(antv4) == 1
        @test_throws ArgumentError set_coordinate_names(antv4, ["u1",  "u2"])
        set_coordinate_names(antv4, ["u"])
        @test_throws ArgumentError set_coordinate_names_of_torus(antv4, ["u"])
        set_coordinate_names_of_torus(antv4, ["u1", "u2"])
        @test is_finalized(antv4) == false
        @test ngens(cox_ring(antv4)) == 1
        @test is_finalized(antv4) == true
    end
    
    @testset "Should fail" begin
        @test_throws ArgumentError ideal_of_linear_relations(antv6)
        @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(antv4)
        @test_throws ArgumentError map_from_torusinvariant_cartier_divisor_group_to_picard_group(antv4)
        @test_throws ArgumentError betti_number(antv4, 0)
        @test_throws ErrorException set_coordinate_names(antv4, ["u"])
        @test_throws ErrorException set_coordinate_names_of_torus(antv4, ["u1",  "u2"])
    end
    
    @testset "Toric ideal" begin
        @test iszero(toric_ideal(antv5)) == true
    end
    
    @testset "Construct from fans" begin
        @test is_trivial(picard_group(antv6)) == true
    end
end