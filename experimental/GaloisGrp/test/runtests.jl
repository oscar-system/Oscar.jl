@testset "SolveByRadicals" begin
  Qx, x = QQ[:x]
  K, r = solve(x^3+3*x+5)
  @test absolute_degree(K) == 12
  @test length(r) == 3

  K, r = solve(x^3+3*x+5; simplify = true)
  @test absolute_degree(K) == 12
  @test length(r) == 3

  K, r = solve(x^3-3)
  @test absolute_degree(K) == 6
  @test length(r) == 3
  @test all(x->x^3 == 3, r)

  K, r = solve(x^3-3; simplify = true)
  @test absolute_degree(K) == 6
  @test length(r) == 3
  @test all(x->x^3 == 3, r)

  Qt, t = rational_function_field(QQ, "t")
  Qtx, x = Qt[:x]
  F, a = function_field(x^6 + 108*t^2 + 108*t + 27)
  s = subfields(F)
  @test length(s) == 4
  G, = galois_group(F)
  @test is_isomorphic(G, symmetric_group(3))
  G,_, k = galois_group(F; overC = true)
  @test is_isomorphic(G, cyclic_group(3))
  @test k isa AbsSimpleNumField && degree(k) == 2

  K, a = cyclotomic_field(3, "a", cached = false)
  G, C = galois_group(K)
  k = fixed_field(C, G)
  @test degree(k) == 1 && k isa AbsSimpleNumField
end

@testset "SubfieldLattice" begin
  Zx, x = ZZ[:x]
  k, a = number_field(swinnerton_dyer(3, x))
  s = subfield_lattice(k)
  @test length(s) == 14
  intersect(s[3], s[4])
  s[3] * s[4]
end
