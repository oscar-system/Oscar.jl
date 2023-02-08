@testset "Hirzebruch surfaces" begin

    F0 = hirzebruch_surface(0)
    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    F5v2 = hirzebruch_surface(5)

    @testset "F0" begin
        @test is_fano(F0) == true
        @test is_projective_space(F0) == false
    end

    @testset "F5" begin
        @test is_normal(F5) == true
        @test is_affine(F5) == false
        @test is_projective(F5) == true
        @test is_smooth(F5) == true
        @test is_complete(F5) == true
        @test has_torusfactor(F5) == false
        @test is_orbifold(F5) == true
        @test is_simplicial(F5) == true
        @test is_gorenstein(F5) == true
        @test is_q_gorenstein(F5) == true
        @test is_fano(F5) == false
        @test dim(F5) == 2
        @test dim_of_torusfactor(F5) == 0
        @test euler_characteristic(F5) == 4
        @test betti_number(F5, 0) == 1
        @test betti_number(F5, 1) == 0
        @test betti_number(F5, 2) == 2
        @test betti_number(F5, 3) == 0
        @test betti_number(F5, 4) == 1
        @test length(affine_open_covering(F5)) == 4
        @test fan(F5).pm_fan.FAN_DIM == 2
        @test rank(torusinvariant_weil_divisor_group(F5)) == 4
        @test rank(character_lattice(F5)) == 2
        @test dim(nef_cone(F5)) == 2
        @test dim(mori_cone(F5)) == 2
        @test rank(domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))) == 2
        @test rank(codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))) == 4
        @test rank(class_group(F5)) == 2
        @test rank(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(F5))) == 2
        @test ngens(cox_ring(F5)) == 4
        @test length(stanley_reisner_ideal(F5).gens) == 2
        @test length(irrelevant_ideal(F5).gens) == 4
        @test is_projective_space(F5) == false
    end

    @testset "Standard constructor" begin
        @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(F5v2)) == matrix(ZZ, [[1, 0], [0, 1], [1, 0], [5, 1]])
        @test transpose(matrix(ZZ,rays(F5v2))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5v2))
        @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5v2)) == character_lattice(F5v2)
        @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5v2)) == torusinvariant_weil_divisor_group(F5v2)
        @test domain(map_from_cartier_divisor_group_to_picard_group(F5v2)) == cartier_divisor_group(F5v2)
        @test codomain(map_from_cartier_divisor_group_to_picard_group(F5v2)) == picard_group(F5v2)
        @test domain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(F5v2)) == cartier_divisor_group(F5v2)
        @test codomain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(F5v2)) == torusinvariant_weil_divisor_group(F5v2)
        @test coordinate_names(F5v2) == ["t1", "x1", "t2", "x2"]
    end
end
