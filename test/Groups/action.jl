@testset "natural stabilizers in permutation groups" begin

  G = symmetric_group(5)
  S = stabilizer(G, 1)
  @test order(S[1]) == 24
  S = stabilizer(G, [1, 2])
  @test order(S[1]) == 6
  S = stabilizer(G, Set([1, 2]))
  @test order(S[1]) == 12

end

@testset "natural stabilizers in matrix groups" begin

  n = 3
  F = GF(2)
  G = general_linear_group(n, F)
  V = AbstractAlgebra.Generic.FreeModule(F, n)
  S = stabilizer(G, gen(V, 1))
  @test order(S[1]) == 24
  S = stabilizer(G, gens(V))
  @test order(S[1]) == 1
  S = stabilizer(G, Set(gens(V)))
  @test order(S[1]) == 6

end
