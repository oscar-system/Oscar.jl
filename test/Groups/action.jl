@testset "natural stabilizers in permutation groups" begin

  G = symmetric_group(5)
  # stabilizer of a vector of integers
  pt = 1
  S = stabilizer(G, pt)
  @test order(S[1]) == 24
  @test S[1] == stabilizer(G, pt, ^)[1]
  # stabilizer of a vector of integers
  pt = [1, 2]
  S = stabilizer(G, pt)
  @test order(S[1]) == 6
  @test S[1] == stabilizer(G, pt, on_tuples)[1]
  # stabilizer of a tuple of integers
  pt = (1, 2)
  S = stabilizer(G, pt)
  @test order(S[1]) == 6
  @test S[1] == stabilizer(G, pt, on_tuples)[1]
  # stabilizer of a set of integers
  pt = Set([1, 2])
  S = stabilizer(G, pt)
  @test order(S[1]) == 12
  @test S[1] == stabilizer(G, pt, on_sets)[1]
  # stabilizer under permutation action on vectors
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

  # larger examples
  G = symmetric_group(100)
  S1, _ = stabilizer(G, [1, 2, 3, 4, 5])
  @test order(S1) == factorial(big(95))
  S2, _ = stabilizer(G, (1, 2, 3, 4, 5))
  @test S2 == S1
  S3, _ = stabilizer(G, Set([1, 2, 3, 4, 5]))
  @test order(S3) == order(S1) * factorial(5)
end

@testset "natural stabilizers in matrix groups" begin
  n = 3
  F = GF(2)
  G = general_linear_group(n, F)
  V = free_module(F, n)
  # stabilizer of a vector
  v = gen(V, 1)
  S = stabilizer(G, v)
  @test order(S[1]) == 24
  @test S[1] == stabilizer(G, v, *)[1]
  @test S[1] == Oscar._stabilizer_generic(G, v, *)[1]
  # stabilizer of a vector of vectors
  S = stabilizer(G, gens(V))
  @test order(S[1]) == 1
  @test S[1] == stabilizer(G, gens(V), on_tuples)[1]
  @test S[1] == Oscar._stabilizer_generic(G, gens(V), on_tuples)[1]
  # stabilizer of a set of vectors
  S = stabilizer(G, Set(gens(V)))
  @test order(S[1]) == 6
  @test S[1] == stabilizer(G, Set(gens(V)), on_sets)[1]
  @test S[1] == Oscar._stabilizer_generic(G, Set(gens(V)), on_sets)[1]
  # stabilizer of a tuple of vectors
  w = (gen(V, 1), gen(V, 2))
  S = stabilizer(G, w)
  @test order(S[1]) == 4
  @test S[1] == stabilizer(G, w, on_tuples)[1]
  @test S[1] == Oscar._stabilizer_generic(G, w, on_tuples)[1]
  # stabilizer of an echelonized matrix
  W, embW = sub(V, [gen(V,1), gen(V,3)])
  m = matrix(embW)
  S = stabilizer(G, m, on_echelon_form_mats)
  @test order(S[1]) == 24
  @test S[1] == Oscar._stabilizer_generic(G, m, on_echelon_form_mats)[1]
  W, embW = sub(V, [])
  m = matrix(embW)
  S = stabilizer(G, m, on_echelon_form_mats)
  @test S[1] == G
  @test S[1] == Oscar._stabilizer_generic(G, m, on_echelon_form_mats)[1]

  # An issue with mismatched OSCAR<->GAP isomorphisms fixed in PR 5150
  E = Oscar.GAPWrap.E
  G_gap = Oscar.GAPWrap.Group(GapObj([
            GapObj([GapObj([1, 0]), GapObj([0, -1])]),
            GapObj(1//2)*GapObj([
              GapObj([ E(24) - E(24)^16 + E(24)^19, E(3)^2 ]),
              GapObj([ -E(3)^2, -E(24) - E(24)^16 - E(24)^19 ])
            ])
          ]))
  G_oscar = Oscar._oscar_group(G_gap)
  K = base_ring(G_oscar)
  V = free_module(K, 2)
  S, _ = stabilizer(G_oscar, [V[1]])
  @test order(S) == 2
end

@testset "action on multivariate polynomials: permutations" begin
  g = symmetric_group(3)
  R, vars = polynomial_ring(QQ, 3);
  (x1, x2, x3) = vars
  f = x1*x2 + x2*x3
  iso = Oscar.iso_oscar_gap(R)
  img = iso(f)

  for p in [cperm(g,1:3), cperm(g,1:2)]
    @test f^p == evaluate(f, permuted(vars, p^-1))
    @test on_indeterminates(img, p) == iso(f^p)
  end

  for x in gens(g), y in gens(g)
    @test on_indeterminates(on_indeterminates(f, x), y) == on_indeterminates(f, x*y)
  end

  I = ideal(R, [x1^2*x2, x2^2])
  orb = orbit(g, on_indeterminates, I)
  @test length(orb) == 6
  II = intersect(collect(orb)...)
  @test length(orbit(g, on_indeterminates, II)) == 1
  @test I^gen(g, 1) == on_indeterminates(I, gen(g, 1))
end

@testset "action on elements of free assoc. algebras: permutations" begin
  g = symmetric_group(3)
  R, vars = free_associative_algebra(QQ, 3);
  (x1, x2, x3) = vars
  f = x1*x2 + x2*x3

  for p in [cperm(g,1:3), cperm(g,1:2)]
    @test f^p == evaluate(f, permuted(vars, p^-1))
  end

  for x in gens(g), y in gens(g)
    @test on_indeterminates(on_indeterminates(f, x), y) == on_indeterminates(f, x*y)
  end
end

@testset "action on multivariate polynomials: matrices" begin
  g = general_linear_group(3, 5)
  R, vars = polynomial_ring(base_ring(g), degree(g))
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

  for x in gens(g), y in gens(g)
    @test on_indeterminates(on_indeterminates(f, x), y) == on_indeterminates(f, x*y)
  end

  I = ideal(R, [x1^2*x2, x2^2])
  orb = orbit(g, on_indeterminates, I)
  @test length(orb) == 186
  II = intersect(collect(orb)...)
  @test length(orbit(g, on_indeterminates, II)) == 1
  @test I^gen(g, 1) == on_indeterminates(I, gen(g, 1))
end

@testset "projective action on lines" begin
  n = 3
  F = GF(5)
  G = general_linear_group(n, F)
  V = free_module(F, n)
  v = gen(V, 1)
  v = on_lines(v, one(G))  # make sure that `v` is normalized
  @test on_lines(2*v, one(G)) == v
  orb = orbit(G, on_lines, v)
  @test length(orb) == 31
  epi = action_homomorphism(orb)
  @test (order(F) - 1) * order(image(epi)[1]) == order(G)
  @test_throws AssertionError on_lines(zero(V), one(G))
end
