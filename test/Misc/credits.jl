# Test that the credited algorithms record their hits in @collect_credits.
@testset "Credits" begin
  @test isempty(@collect_credits 1 + 1)

  F3 = GF(3)
  N1 = matrix(F3, 3, 3, [0, 1, 0, 2, 0, 0, 0, 0, 2])
  N2 = matrix(F3, 3, 3, [2, 0, 0, 0, 2, 0, 0, 0, 2])
  RG = invariant_ring(N1, N2)
  cred = @collect_credits primary_invariants(RG)
  @test only(keys(cred)) == Oscar.CREDIT_Kem99
  @test occursin("@Article{Kem99,", AbstractAlgebra.Credits.bibtex(cred))
  @test isempty(@collect_credits primary_invariants(RG))   # cached

  cred = @collect_credits begin
    collect(partitions(7, 3))
    collect(partitions(100, 7, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))
  end
  @test only(keys(cred)) == Oscar.CREDIT_RJ76
  @test length(cred[Oscar.CREDIT_RJ76]) == 2

  cred = @collect_credits begin
    pk, sk = lwe_keygen(lwe_toy)
    lwe_decrypt(lwe_toy, sk, lwe_encrypt(lwe_toy, pk, UInt8[0x2a]))
  end
  @test only(keys(cred)) == Oscar.CREDIT_Reg05
  @test sum(values(cred[Oscar.CREDIT_Reg05])) == 3
end
