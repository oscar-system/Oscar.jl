@testset "Line bundles" begin
  dP1 = del_pezzo_surface(NormalToricVariety, 1)
  dP3 = del_pezzo_surface(NormalToricVariety, 3)

  l = toric_line_bundle(dP3, [1, 2, 3, 4])
  l2 = canonical_bundle(dP3)
  l3 = anticanonical_bundle(dP3)
  l4 = toric_line_bundle(dP3, trivial_divisor(dP3))
  l5 = toric_line_bundle(dP1, [ZZRingElem(1), ZZRingElem(2)])

  @testset "Should fail" begin
    @test_throws ArgumentError l * l5
  end

  @testset "Basic properties" begin
    @test is_trivial(l) == false
    @test is_basepoint_free(l) == false
    @test is_ample(l) == false
    @test is_very_ample(l) == false
  end

  @testset "Basic attributes" begin
    @test degree(l) == -6
    @test degree(l^(-1)) == 6
    @test degree(l*l) == -12
    @test picard_class(l).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
    @test dim(toric_variety(l)) == 2
  end

  @testset "Arithmetic" begin
    @test (l == l5) == false
    @test (l == l2) == false
    @test (l2 * l3 == structure_sheaf(dP3)) == true
    @test (l * l4 * inv(l) == structure_sheaf(dP3)) == true
  end
end
