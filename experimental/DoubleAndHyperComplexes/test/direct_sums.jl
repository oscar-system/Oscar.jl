@testset "direct sums of complexes" begin
  R, (x, y) = QQ[:x, :y]
  F = FreeMod(R, 1)
  I, _ = ideal(gens(R))*F
  res, _ = free_resolution(Oscar.SimpleFreeResolution, I)
  C = Oscar.ZeroDimensionalComplex(I)
  D = direct_sum([res, res])
  @test ngens(D[0]) == 4
  @test ngens(D[1]) == 2
  E = direct_sum([C, C])
  @test ngens(E[()]) == 4
  hC = hom(res, C)
  # The following combination of input does not produce a type for the vector 
  # from which the chain- and morphism-type of the direct sum can still be inferred
  @test_throws ErrorException direct_sum([res, hC])
  # In such cases the user has to specify the types explicitly to the internal constructor. 
  DD = Oscar.DirectSumComplex{OFPModule, OFPModuleHom}([res, hC])
  @test ngens(DD[0]) == 6
  @test !is_zero(DD[-1])
  @test !is_zero(DD[1])
  DD = Oscar.DirectSumComplex{OFPModule, OFPModuleHom}([res[2:3], hC])
  @test !Oscar.can_compute_index(DD, 7) # none of the summands is defined there
end

