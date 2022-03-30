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

  # a more complex example
  G = symmetric_group(14)
  gens = [ cperm([1,10], [2,12,13,7,8,14], [3,4,9,6,5,11]),
           cperm([1,11,6,13,12,10,2,8,4,3]) ]
  H = sub(G, gens)[1]

  # stabilizer should work
  K = stabilizer(H, 1)[1]
  @test order(K) == 46080

  # bugfix test: stabilizer should still work after group size is known
  @test order(H) == 645120
  @test K == stabilizer(H, 1)[1]

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
