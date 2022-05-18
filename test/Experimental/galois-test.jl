@testset "Experimental.galois_group" begin

  Zx, x = ZZ["x"]
  k, a = number_field(x^5-2)
  G, C = galois_group(k)
  @test transitive_group_identification(G) == (5, 3)

  U = trivial_subgroup(G)[1]
  L = fixed_field(C, U)
  @test degree(L) == order(G)
  @test length(roots(k.pol, L)) == 5

  R, x = PolynomialRing(QQ, "x")
  pol = x^6 - 366*x^4 - 878*x^3 + 4329*x^2 + 14874*x + 10471
  g, C = galois_group(pol)
  @test order(g) == 18

  s = symmetric_group(4)
  g = sub(s, [s([2,1,4,3]), s([3,4,1,2])])[1]
  act = Oscar.GaloisGrp.action_on_blocks(g, [1, 2])
  @test order(image(act)[1]) == 2

end
