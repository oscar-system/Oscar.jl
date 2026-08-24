# Test defining data, arithmetic, and representative-independent degree of toric line bundles.
@testset "Line bundles" begin
  dP1 = del_pezzo_surface(NormalToricVariety, 1)
  dP3 = del_pezzo_surface(NormalToricVariety, 3)

  l = toric_line_bundle(dP3, [1, 2, 3, 4])
  l2 = canonical_bundle(dP3)
  l3 = anticanonical_bundle(dP3)
  l4 = toric_line_bundle(dP3, trivial_divisor(dP3))
  l5 = toric_line_bundle(dP1, [ZZRingElem(1), ZZRingElem(2)])
  P2 = projective_space(NormalToricVariety, 2)
  l6 = toric_line_bundle(P2, [2])

  @testset "Should fail due to bad arguments (toric line bundles)" begin
    @test_throws ArgumentError l * l5
  end

  @testset "Basic properties" begin
    @test is_trivial(l) == false
    @test is_basepoint_free(l) == false
    @test is_ample(l) == false
    @test is_very_ample(l) == false
  end

  @testset "Basic attributes" begin
    @test_throws ArgumentError degree(l)
    @test degree(l6) == 2
    @test degree(l6^(-1)) == -2
    @test degree(l6 * l6) == 4
    @test picard_class(l).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
    @test dim(toric_variety(l)) == 2
  end

  @testset "Degree is independent of the divisor representative" begin
    WPS = weighted_projective_space(NormalToricVariety, [2, 3, 1])
    principal_divisor = divisor_of_character(WPS, [1, 2])
    principal_bundle = toric_line_bundle(principal_divisor)
    trivial_bundle = trivial_line_bundle(WPS)
    @test principal_bundle == trivial_bundle
    @test degree(principal_bundle) == degree(trivial_bundle) == 0
  end

  @testset "Arithmetic" begin
    @test (l == l5) == false
    @test (l == l2) == false
    @test (l2 * l3 == structure_sheaf(dP3)) == true
    @test (l * l4 * inv(l) == structure_sheaf(dP3)) == true
  end
end
