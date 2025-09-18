@testset "Projectivization of direct sum of line bundles on toric space" begin
  P1 = projective_space(NormalToricVariety, 1)
  @testset "Test of the Hirzebruch surfaces from 0 to 10" begin
    for a in 0:10
      Da = toric_divisor(P1, [a, 0])
      D0 = toric_divisor(P1, [0, 0])
      X = projectivization(D0, Da)
      @test is_smooth(X) == true
      @test (a < 2) == is_fano(X)
      @test torsion_free_rank(picard_group_with_map(X)[1]) == 2
      @test integrate(cohomology_class(anticanonical_divisor(X))^dim(X)) == integrate(
        cohomology_class(anticanonical_divisor(hirzebruch_surface(NormalToricVariety, a)))^2
      )
    end
  end

  # We will test the rank of the Picard group and the constant (-K_X)^3
  # Our reference is https://www.fanography.info/toric

  # Let us start from the projective bundle over P2
  P2 = projective_space(NormalToricVariety, 2)
  D0 = toric_divisor(P2, [0, 0, 0])
  X1_P2 = projectivization(D0, D0)
  D1 = toric_divisor(P2, [0, 0, 1])
  X2_P2 = projectivization(D0, D1)
  D2 = toric_divisor(P2, [0, 2, 0])
  X3_P2 = projectivization(D0, D2)
  @testset "Test of some Fano projective bundles of dimension 3 over P2" begin
    @test torsion_free_rank(picard_group_with_map(X1_P2)[1]) == 2
    @test is_fano(X1_P2) == true
    @test integrate(cohomology_class(anticanonical_divisor(X1_P2))^dim(X1_P2)) == 54
    @test torsion_free_rank(picard_group_with_map(X2_P2)[1]) == 2
    @test is_fano(X2_P2) == true
    @test integrate(cohomology_class(anticanonical_divisor(X2_P2))^dim(X2_P2)) == 56
    @test torsion_free_rank(picard_group_with_map(X3_P2)[1]) == 2
    @test is_fano(X3_P2) == true
    @test integrate(cohomology_class(anticanonical_divisor(X3_P2))^dim(X3_P2)) == 62
  end

  # Projective bundle over F1 = hirzebruch_surface(1)
  D0 = toric_divisor(P1, [0, 0])
  D1 = toric_divisor(P1, [-1, 0])
  F1 = projectivization(D0, D1)
  X1_F1 = projectivization(toric_divisor(F1, [1, 1, 0, 0]), toric_divisor(F1, [1, 1, 0, 0]))
  l1_F1 = toric_line_bundle(F1, [0, 1])
  X2_F1 = projectivization(toric_line_bundle(F1, [0, 0]), l1_F1)
  @testset "Test of some Fano projective bundles of dimension 3 over F1" begin
    @test torsion_free_rank(picard_group_with_map(X1_F1)[1]) == 3
    @test is_fano(X1_F1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X1_F1))^dim(X1_F1)) == 48
    @test integrate(cohomology_class(l1_F1)^2) == 1
    @test torsion_free_rank(picard_group_with_map(X2_F1)[1]) == 3
    @test is_fano(X2_F1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X2_F1))^dim(X2_F1)) == 50
  end

  P1xP1 = projectivization(D0, D0)
  O00 = toric_line_bundle(P1xP1, [0, 0])
  O10 = toric_line_bundle(P1xP1, [1, 0])
  O01 = toric_line_bundle(P1xP1, [0, 1])
  O11 = toric_line_bundle(P1xP1, [1, 1])
  X1_P1xP1 = projectivization(O10, O01)
  X2_P1xP1 = projectivization(O10, O10)
  X3_P1xP1 = projectivization(O00, O11)
  @testset "Test of some Fano projective bundles of dimension 3 over P1 * P1" begin
    @test torsion_free_rank(picard_group_with_map(X1_P1xP1)[1]) == 3
    @test is_fano(X1_P1xP1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X1_P1xP1))^dim(X1_P1xP1)) == 44
    @test torsion_free_rank(picard_group_with_map(X2_P1xP1)[1]) == 3
    @test is_fano(X2_P1xP1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X2_P1xP1))^dim(X2_P1xP1)) == 48
    @test torsion_free_rank(picard_group_with_map(X3_P1xP1)[1]) == 3
    @test is_fano(X3_P1xP1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X3_P1xP1))^dim(X3_P1xP1)) == 52
  end
end
