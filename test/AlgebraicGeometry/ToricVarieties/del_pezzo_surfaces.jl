@testset "del Pezzo surfaces" begin

  dP0 = del_pezzo_surface(NormalToricVariety, 0)
  dP1 = del_pezzo_surface(NormalToricVariety, 1)
  dP2 = del_pezzo_surface(NormalToricVariety, 2)
  dP3 = del_pezzo_surface(NormalToricVariety, 3)

  @testset "Should fail" begin
    @test_throws ArgumentError del_pezzo_surface(NormalToricVariety, -1)
    @test_throws ArgumentError del_pezzo_surface(NormalToricVariety, 4)
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
    @test torsion_free_rank(torusinvariant_cartier_divisor_group(dP1)) == 4
    @test torsion_free_rank(picard_group(dP1)) == 2
    @test transpose(matrix(ZZ,rays(dP1))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP1))
    @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP1)) == character_lattice(dP1)
    @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP1)) == torusinvariant_weil_divisor_group(dP1)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP1)) == torusinvariant_cartier_divisor_group(dP1)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP1)) == picard_group(dP1)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(dP1)) == torusinvariant_cartier_divisor_group(dP1)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(dP1)) == torusinvariant_weil_divisor_group(dP1)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_class_group(dP1)) == torusinvariant_cartier_divisor_group(dP1)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_class_group(dP1)) == class_group(dP1)
    @test domain(map_from_picard_group_to_class_group(dP1)) == picard_group(dP1)
    @test codomain(map_from_picard_group_to_class_group(dP1)) == class_group(dP1)
    @test gorenstein_index(dP1) == 1
    @test picard_index(dP1) == 1
    @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(dP1)) == matrix(ZZ, [[1, 1], [1, 1], [1, 0], [0, -1]])
    @test coordinate_names(dP1) == ["x1", "x2", "x3", "e1"]
  end
  
  @testset "Basic attributes of dP2" begin
    @test length(torusinvariant_prime_divisors(dP2)) == 5
    @test torsion_free_rank(torusinvariant_cartier_divisor_group(dP2)) == 5
    @test torsion_free_rank(picard_group(dP2)) == 3
    @test transpose(matrix(ZZ,rays(dP2))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP2))
    @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP2)) == character_lattice(dP2)
    @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP2)) == torusinvariant_weil_divisor_group(dP2)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP2)) == torusinvariant_cartier_divisor_group(dP2)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP2)) == picard_group(dP2)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(dP2)) == torusinvariant_cartier_divisor_group(dP2)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(dP2)) == torusinvariant_weil_divisor_group(dP2)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_class_group(dP2)) == torusinvariant_cartier_divisor_group(dP2)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_class_group(dP2)) == class_group(dP2)
    @test domain(map_from_picard_group_to_class_group(dP2)) == picard_group(dP2)
    @test codomain(map_from_picard_group_to_class_group(dP2)) == class_group(dP2)
    @test gorenstein_index(dP2) == 1
    @test picard_index(dP2) == 1
    @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(dP2)) == matrix(ZZ, [[1, 1, 1], [1, 1, 0], [1, 0, 1], [0, -1, 0], [0, 0, -1]])
    @test coordinate_names(dP2) == ["x1", "x2", "x3", "e1", "e2"]
  end
  
  @testset "Basic attributes of dP3" begin
    @test length(torusinvariant_prime_divisors(dP3)) == 6
    @test torsion_free_rank(torusinvariant_cartier_divisor_group(dP3)) == 6
    @test torsion_free_rank(picard_group(dP3)) == 4
    @test transpose(matrix(ZZ,rays(dP3))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP3))
    @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP3)) == character_lattice(dP3)
    @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(dP3)) == torusinvariant_weil_divisor_group(dP3)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP3)) == torusinvariant_cartier_divisor_group(dP3)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP3)) == picard_group(dP3)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(dP3)) == torusinvariant_cartier_divisor_group(dP3)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(dP3)) == torusinvariant_weil_divisor_group(dP3)
    @test domain(map_from_torusinvariant_cartier_divisor_group_to_class_group(dP3)) == torusinvariant_cartier_divisor_group(dP3)
    @test codomain(map_from_torusinvariant_cartier_divisor_group_to_class_group(dP3)) == class_group(dP3)
    @test domain(map_from_picard_group_to_class_group(dP3)) == picard_group(dP3)
    @test codomain(map_from_picard_group_to_class_group(dP3)) == class_group(dP3)
    @test gorenstein_index(dP3) == 1
    @test picard_index(dP3) == 1
    @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(dP3)) == matrix(ZZ, [[1, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    @test coordinate_names(dP3) == ["x1", "x2", "x3", "e1", "e2", "e3"]
  end

  @testset "nef cone" begin
    nc = nef_cone(dP3)
    rnc = matrix(ZZ, rays(nc))
    lnc = [toric_line_bundle(dP3, vec([ZZ(k) for k in rnc[l,:]])) for l in 1:nrows(rnc)]
    cnc = [cohomology_class(l) for l in lnc]
    dnc = [toric_divisor(l) for l in lnc]
    for d in dnc
      @test is_nef(d)
    end
    @test is_nef(dnc[1] + dnc[2])
  end
end
