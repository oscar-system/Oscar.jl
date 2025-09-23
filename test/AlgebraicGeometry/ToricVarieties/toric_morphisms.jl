@testset "Toric morphisms" begin
  source = projective_space(NormalToricVariety, 1)
  target = hirzebruch_surface(NormalToricVariety, 2)
  tm1 = toric_morphism(source, matrix(ZZ, [[0, 1]]), target)
  tm2 = toric_identity_morphism(target)
  tm3 = toric_identity_morphism(hirzebruch_surface(NormalToricVariety, 4))
  p1p1 = normal_toric_variety(cube(2))
  p2 = projective_space(NormalToricVariety, 2)

  @testset "Argument errors for toric morphisms" begin
    @test_throws ArgumentError tm1 * tm1
    @test_throws ArgumentError toric_morphism(p1p1, [1 0; 0 1], p2)
    @test_throws ArgumentError toric_morphism(p2, [[1, 0], [0, 1]], p1p1)
  end

  @testset "Arithmetic of toric morphisms" begin
    @test (tm1 == tm2) == false
    @test (tm2 + tm2 == 2*tm2) == true
    @test (tm2 + tm2 == ZZRingElem(2)*tm2) == true
    @test (tm2 * tm2 == tm2) == true
  end

  @testset "Basic attributes of toric morphisms" begin
    @test codomain(tm1) === domain(tm2)
    @test matrix(lattice_homomorphism(tm2)) == matrix(ZZ, [[1, 0], [0, 1]])
    @test morphism_on_torusinvariant_weil_divisor_group(tm2).map == identity_matrix(ZZ, 4)
    @test morphism_on_torusinvariant_cartier_divisor_group(tm2).map ==
      identity_matrix(ZZ, 4)
    @test lattice_homomorphism(morphism_from_cox_variety(source)).map ==
      matrix(ZZ, [[1], [-1]])
    @test is_affine(cox_variety(source)) == false
    @test matrix(morphism_on_class_group(tm3)) == matrix(ZZ, [[1, 0], [0, 1]])
    @test matrix(morphism_on_picard_group(tm3)) == matrix(ZZ, [[1, 0], [0, 1]])
  end
end
