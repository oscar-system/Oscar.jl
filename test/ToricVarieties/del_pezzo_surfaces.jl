using Oscar
using Test

@testset "del Pezzo surfaces (set_attributes = $set_attributes)" for set_attributes in [true, false]

    dP0 = del_pezzo_surface(0; set_attributes)
    dP1 = del_pezzo_surface(1; set_attributes)
    dP2 = del_pezzo_surface(2; set_attributes)
    dP3 = del_pezzo_surface(3; set_attributes)

    @testset "Should fail" begin
        @test_throws ArgumentError del_pezzo_surface(-1; set_attributes)
        @test_throws ArgumentError del_pezzo_surface(4; set_attributes)
    end

    @testset "Basic properties" begin
        @test is_projective_space(dP0) == true
        @test is_projective_space(dP1) == false
        @test is_projective_space(dP2) == false
        @test is_projective_space(dP3) == false
    end
    
    @testset "Basic attributes of dP0" begin
        @test length(torusinvariant_prime_divisors(dP0)) == 3
    end
    
    @testset "Basic attributes of dP1" begin
        @test length(torusinvariant_prime_divisors(dP1)) == 4
        @test rank(torusinvariant_cartier_divisor_group(dP1)) == 4
        @test rank(picard_group(dP1)) == 2
        @test transpose(matrix(ZZ,rays(dP1))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP1))
        @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP1)) == character_lattice(dP1)
        @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP1)) == torusinvariant_weil_divisor_group(dP1)
        @test domain(map_from_cartier_divisor_group_to_picard_group(dP1)) == cartier_divisor_group(dP1)
        @test codomain(map_from_cartier_divisor_group_to_picard_group(dP1)) == picard_group(dP1)
        @test domain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(dP1)) == cartier_divisor_group(dP1)
        @test codomain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(dP1)) == torusinvariant_weil_divisor_group(dP1)
        if set_attributes
            @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(dP1)) == matrix(ZZ, [[1, 1], [1, 1], [1, 0], [0, -1]])
            @test coordinate_names(dP1) == ["x1", "x2", "x3", "e1"]
        end
    end
    
    @testset "Basic attributes of dP2" begin
        @test length(torusinvariant_prime_divisors(dP2)) == 5
        @test rank(torusinvariant_cartier_divisor_group(dP2)) == 5
        @test rank(picard_group(dP2)) == 3
        @test transpose(matrix(ZZ,rays(dP2))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP2))
        @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP2)) == character_lattice(dP2)
        @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP2)) == torusinvariant_weil_divisor_group(dP2)
        @test domain(map_from_cartier_divisor_group_to_picard_group(dP2)) == cartier_divisor_group(dP2)
        @test codomain(map_from_cartier_divisor_group_to_picard_group(dP2)) == picard_group(dP2)
        @test domain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(dP2)) == cartier_divisor_group(dP2)
        @test codomain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(dP2)) == torusinvariant_weil_divisor_group(dP2)
        if set_attributes
            @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(dP2)) == matrix(ZZ, [[1, 1, 1], [1, 1, 0], [1, 0, 1], [0, -1, 0], [0, 0, -1]])
            @test coordinate_names(dP2) == ["x1", "x2", "x3", "e1", "e2"]
        end
    end
    
    @testset "Basic attributes of dP3" begin
        @test length(torusinvariant_prime_divisors(dP3)) == 6
        @test rank(torusinvariant_cartier_divisor_group(dP3)) == 6
        @test rank(picard_group(dP3)) == 4
        @test transpose(matrix(ZZ,rays(dP3))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP3))
        @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP3)) == character_lattice(dP3)
        @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP3)) == torusinvariant_weil_divisor_group(dP3)
        @test domain(map_from_cartier_divisor_group_to_picard_group(dP3)) == cartier_divisor_group(dP3)
        @test codomain(map_from_cartier_divisor_group_to_picard_group(dP3)) == picard_group(dP3)
        @test domain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(dP3)) == cartier_divisor_group(dP3)
        @test codomain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(dP3)) == torusinvariant_weil_divisor_group(dP3)
        if set_attributes
            @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(dP3)) == matrix(ZZ, [[1, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
            @test coordinate_names(dP3) == ["x1", "x2", "x3", "e1", "e2", "e3"]
        end
    end
end
