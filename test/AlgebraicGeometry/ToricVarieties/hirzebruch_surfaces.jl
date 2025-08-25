@testset "Hirzebruch surfaces" begin
  
  F0 = hirzebruch_surface(NormalToricVariety, 0)
  F5 = hirzebruch_surface(NormalToricVariety, 5)
  
  @testset "F0" begin
    @test is_fano(F0) == true
    @test is_projective_space(F0) == false
  end
  
  @testset "Properties of F5" begin
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
    @test is_projective_space(F5) == false
  end
  
  @testset "Attributes of F5" begin
    @test dim(F5) == 2
    @test dim_of_torusfactor(F5) == 0
    @test euler_characteristic(F5) == 4
    @test betti_number(F5, 0) == 1
    @test betti_number(F5, 1) == 0
    @test betti_number(F5, 2) == 2
    @test betti_number(F5, 3) == 0
    @test betti_number(F5, 4) == 1
    @test length(affine_open_covering(F5)) == 4
    @test dim(polyhedral_fan(F5)) == 2
    @test torsion_free_rank(torusinvariant_weil_divisor_group(F5)) == 4
    @test torsion_free_rank(character_lattice(F5)) == 2
    @test ngens(cox_ring(F5)) == 4
    @test length(stanley_reisner_ideal(F5).gens) == 2
    @test length(irrelevant_ideal(F5).gens) == 4
    @test dim(nef_cone(F5)) == 2
    @test dim(mori_cone(F5)) == 2
    @test torsion_free_rank(domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))) == 2
    @test torsion_free_rank(codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))) == 4
    @test torsion_free_rank(class_group_with_map(F5)[1]) == 2
    @test torsion_free_rank(codomain(map_from_torusinvariant_weil_divisor_group_to_class_group(F5))) == 2
    @test transpose(matrix(ZZ,rays(F5))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5))
    @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5)) == character_lattice(F5)
    @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(F5)) == torusinvariant_weil_divisor_group(F5)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(F5)) == torusinvariant_cartier_divisor_group(F5)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(F5)) == picard_group_with_map(F5)[1]
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(F5)) == torusinvariant_cartier_divisor_group(F5)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(F5)) == torusinvariant_weil_divisor_group(F5)
    @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(F5)) == matrix(ZZ, [[1, 0], [0, 1], [1, 0], [5, 1]])
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_class_group(F5)) == torusinvariant_cartier_divisor_group(F5)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_class_group(F5)) == class_group_with_map(F5)[1]
    @test domain(map_from_picard_group_to_class_group(F5)) == picard_group_with_map(F5)[1]
    @test codomain(map_from_picard_group_to_class_group(F5)) == class_group_with_map(F5)[1]
    @test gorenstein_index(F5) == 1
    @test picard_index(F5) == 1
    @test coordinate_names(F5) == ["t1", "x1", "t2", "x2"]
  end
end
