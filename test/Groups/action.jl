@testset "natural stabilizers in permutation groups" begin

  G = symmetric_group(5)
  S = stabilizer(G, 1)
  @test order(S[1]) == 24
  @test S[1] == stabilizer(G, 1, ^)[1]
  pt = [1, 2]
  S = stabilizer(G, pt)
  @test order(S[1]) == 6
  @test S[1] == stabilizer(G, pt, on_tuples)[1]
  S = stabilizer(G, Set(pt))
  @test order(S[1]) == 12
  @test S[1] == stabilizer(G, Set(pt), on_sets)[1]
  l = [1, 1, 2, 2, 3]
  S = stabilizer(G, l, permuted)
  @test order(S[1]) == 4

end

@testset "natural stabilizers in matrix groups" begin

  n = 3
  F = GF(2)
  G = general_linear_group(n, F)
  V = AbstractAlgebra.Generic.FreeModule(F, n)
  v = gen(V, 1)
  S = stabilizer(G, v)
  @test order(S[1]) == 24
  @test S[1] == stabilizer(G, v, *)[1]
  S = stabilizer(G, gens(V))
  @test order(S[1]) == 1
  @test S[1] == stabilizer(G, gens(V), on_tuples)[1]
  S = stabilizer(G, Set(gens(V)))
  @test order(S[1]) == 6
  @test S[1] == stabilizer(G, Set(gens(V)), on_sets)[1]

end
