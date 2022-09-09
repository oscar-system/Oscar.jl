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

@testset "action on multivariate polynomials: permutations" begin
  g = symmetric_group(3)
  R, vars = PolynomialRing(QQ, 3);
  (x1, x2, x3) = vars
  f = x1*x2 + x2*x3
  iso = Oscar.iso_oscar_gap(R)
  img = iso(f)

  for p in [g(cperm(1:3)), g(cperm(1:2))]
    @test f^p == evaluate(f, permuted(vars, p^-1))
    @test on_indeterminates(img, p) == iso(f^p)
  end

  for x in gens(g)
    for y in gens(g)
      @test on_indeterminates(on_indeterminates(f, x), y) == on_indeterminates(f, x*y)
    end
  end

  I = ideal(R, [x1^2*x2, x2^2])
  orb = orbit(g, on_indeterminates, I)
  @test length(orb) == 6
  II = intersect(collect(orb)...)
  @test length(orbit(g, on_indeterminates, II)) == 1
  @test I^gen(g, 1) == on_indeterminates(I, gen(g, 1))
end

@testset "action on multivariate polynomials: matrices" begin
  g = general_linear_group(3, 5)
  R, vars = PolynomialRing(base_ring(g), degree(g))
  (x1, x2, x3) = vars
  f = x1*x2 + x2*x3

  # permutation matrix
  p = cperm(1:3)
  m = g(permutation_matrix(base_ring(g), p))
  @test f^p == x1*x3 + x2*x3
  @test f^m == f^p
  iso = Oscar.iso_oscar_gap(R)
  img = iso(f)
  @test on_indeterminates(img, m) == iso(f^m)

  # non-permutation matrix
  m = g(matrix(base_ring(g), 3, 3, [3, 0, 2, 4, 0, 0, 0, 4, 0]))
  @test f^m == 2*x1^2 + x1*x2 + 3*x1*x3

  for x in gens(g)
    for y in gens(g)
      @test on_indeterminates(on_indeterminates(f, x), y) == on_indeterminates(f, x*y)
    end
  end

  I = ideal(R, [x1^2*x2, x2^2])
  orb = orbit(g, on_indeterminates, I)
  @test length(orb) == 186
  II = intersect(collect(orb)...)
  @test length(orbit(g, on_indeterminates, II)) == 1
  @test I^gen(g, 1) == on_indeterminates(I, gen(g, 1))
end
