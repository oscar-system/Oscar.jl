@testset "Presentation as affine algebra" begin
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 0, 0, 1, 1, 0, 0 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])

  RG = invariant_ring(M1, M2)
  A, AtoR = affine_algebra(RG)
  @test [ AtoR(x) for x in gens(A) ] == fundamental_invariants(RG)
  @test is_injective(AtoR)

  RG = invariant_ring(M1, M2)
  A, AtoR = affine_algebra(RG, algo_rels = :linear_algebra)
  @test [ AtoR(x) for x in gens(A) ] == fundamental_invariants(RG)
  @test is_injective(AtoR)
end
