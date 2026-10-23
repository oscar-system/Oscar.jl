# Test class accessors, Cartier conversion, and arithmetic compatibility checks.
@testset "Torus-invariant divisor classes" begin
  F5 = hirzebruch_surface(NormalToricVariety, 5)
  dP3 = del_pezzo_surface(NormalToricVariety, 3)
  P2 = projective_space(NormalToricVariety, 2)

  rs = [
    [0, 1, 0],
    [-1, 1, 1],
    [-1, 0, 0],
    [-1, -2, -2],
    [0, -1, -1],
    [0, 0, 1],
    [-1, -1, -1],
    [1, 0, 0],
  ]
  cs = incidence_matrix([
    [1, 2, 3],
    [1, 4, 5],
    [3, 6, 7],
    [4, 5, 6],
    [1, 4, 7],
    [5, 6, 8],
    [1, 6, 8],
    [2, 3, 6],
    [1, 5, 8],
    [1, 3, 7],
    [1, 2, 6],
    [4, 6, 7],
  ])
  test_space = normal_toric_variety(cs, rs)

  DC = toric_divisor_class(F5, [ZZRingElem(0), ZZRingElem(0)])
  DC2 = toric_divisor_class(F5, [1, 2])
  DC3 = toric_divisor_class(dP3, [4, 3, 2, 1])
  DC4 = canonical_divisor_class(dP3)
  DC5 = anticanonical_divisor_class(dP3)
  DC6 = trivial_divisor_class(dP3)
  DC7 = toric_divisor_class(P2, [1])
  DC8 = toric_divisor_class(P2, [-1])

  @testset "Basic properties" begin
    @test is_trivial(toric_divisor(DC2)) == false
    @test is_effective(DC7) == true
    @test is_effective(DC8) == false
    @test is_effective(anticanonical_divisor_class(test_space)) == true
  end

  @testset "Basic attributes" begin
    @test torsion_free_rank(parent(divisor_class(DC2))) == 2
    @test coefficients(DC2) == coefficients(toric_divisor(DC2))
    @test picard_class(DC7) == picard_class(toric_line_bundle(DC7))
    @test map_from_picard_group_to_class_group(P2)(picard_class(DC7)) == divisor_class(DC7)
    @test dim(toric_variety(DC2)) == 2
  end

  @testset "Arithmetic" begin
    @test is_trivial(ZZRingElem(2) * DC + DC2) == false
    @test is_trivial(2 * DC - DC2) == false
    @test (DC == DC2) == false
    @test (DC4 - DC5 == DC6) == false
    @test (DC == DC3) == false
  end

  @testset "Should fail due to incompatible or non-Cartier classes" begin
    another_F5 = hirzebruch_surface(NormalToricVariety, 5)
    another_DC = toric_divisor_class(another_F5, [1, 2])
    @test_throws ArgumentError DC2 + another_DC
    @test_throws ArgumentError DC2 - another_DC
    WPS = weighted_projective_space(NormalToricVariety, [2, 3, 1])
    non_cartier_class = toric_divisor_class(toric_divisor(WPS, [1, 0, 0]))
    @test_throws ArgumentError picard_class(non_cartier_class)
    @test_throws ArgumentError toric_line_bundle(non_cartier_class)
  end
end
