@testset "Experimental.galois_group" begin

  Zx, x = ZZ["x"]
  k, a = number_field(x^5-2)
  G, C = galois_group(k)
  @test transitive_identification(G) == 3

  U = trivial_subgroup(G)[1]
  L = fixed_field(C, U)
  @test degree(L) == order(G)
  @test length(roots(k.pol, L)) == 5
end
