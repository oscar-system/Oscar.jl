@testset "Invariant Theory of SL_m" begin
  #h is invariant if h(M*x) - h(x) vanishes modulo the defining ideal of the group
  function is_invariant(h::MPolyRingElem, r::Oscar.RepresentationLinearlyReductiveGroup)
    M = representation_matrix(r)
    S = base_ring(M)
    Sx, x = polynomial_ring(S, :x => 1:ncols(M); cached=false)
    hx = map_coefficients(S, forget_grading(h); parent=Sx)
    d = evaluate(hx, M * x) - hx
    return all(c -> is_zero(normal_form(c, defining_ideal(group(r)))), coefficients(d))
  end

  S, z = polynomial_ring(QQ, :z => (1:2, 1:2))
  G = linearly_reductive_group(:SL, 2, S)
  @test Oscar.group_type(G) == :SL
  @test Oscar.group_dim(G) == 2
  @test defining_ideal(G) == ideal([z[1, 1] * z[2, 2] - z[2, 1] * z[1, 2] - 1])
  @test Oscar.natural_representation(G) == matrix(S, z)
  @test_throws ArgumentError linearly_reductive_group(:GL, 2, QQ)
  @test_throws ArgumentError linearly_reductive_group(:SL, 2, GF(5))

  rep1 = representation_on_forms(G, 2)
  @test ncols(representation_matrix(rep1)) == 3
  @test group(rep1) == G
  @test vector_space_dim(rep1) == 3
  @test degree(rep1) == 3
  @test sprint(show, rep1) == "Representation of degree 3 of linearly reductive group"

  R_rep1 = invariant_ring(rep1)
  FI_rep1 = fundamental_invariants(R_rep1)
  @test length(FI_rep1) == 1

  #the null cone ideal is Hilbert's ideal, in the polynomial ring of the invariant ring
  N = null_cone_ideal(R_rep1)
  @test base_ring(N) === polynomial_ring(R_rep1)
  @test N == ideal(FI_rep1)

  #same rep mat as in rep1, only without multinomial coefficients.
  M = matrix(
    S,
    3,
    3,
    [
      z[1, 1]^2 z[1, 1]*z[2, 1] z[2, 1]^2;
      2*z[1, 1]*z[1, 2] z[1, 1] * z[2, 2]+z[2, 1] * z[1, 2] 2*z[2, 1]*z[2, 2];
      z[1, 2]^2 z[1, 2]*z[2, 2] z[2, 2]^2
    ],
  )
  rep2 = representation_reductive_group(G, M)
  R_rep2 = invariant_ring(R_rep1.poly_ring, rep2)
  FI_rep2 = fundamental_invariants(R_rep2)
  x = gens(R_rep1.poly_ring)
  @test FI_rep2 == [-4 * x[1] * x[3] + x[2]^2]

  #Reynolds operator on inhomogeneous polynomials
  h = FI_rep1[1]
  @test reynolds_operator(R_rep1, h^2 + h + 1) == h^2 + h + 1
  @test is_zero(reynolds_operator(R_rep1, zero(h)))

  #representation matrix whose rows are not homogeneous
  M = block_diagonal_matrix([
    representation_matrix(representation_on_forms(G, 1)), representation_matrix(rep1)
  ])
  P = identity_matrix(S, 5)
  P[1, 3] = 1
  rep_mixed = representation_reductive_group(G, inv(P) * M * P)
  FI_mixed = fundamental_invariants(invariant_ring(rep_mixed))
  @test sort(map(total_degree, FI_mixed)) == [2, 3]
  @test all(h -> is_invariant(h, rep_mixed), FI_mixed)

  #fields other than QQ
  Qt, _ = rational_function_field(QQ, :t)
  for K in [quadratic_field(2)[1], Qt, algebraic_closure(QQ)]
    GK = linearly_reductive_group(:SL, 2, K)
    repK = representation_on_forms(GK, 2)
    FI_K = fundamental_invariants(invariant_ring(repK))
    @test length(FI_K) == 1
    @test is_invariant(FI_K[1], repK)
  end

  #no invariants of positive degree
  R_V1 = invariant_ring(representation_on_forms(G, 1))
  @test is_empty(fundamental_invariants(R_V1))
  @test is_zero(null_cone_ideal(R_V1))
  @test base_ring(null_cone_ideal(R_V1)) === polynomial_ring(R_V1)

  #the representation matrix maps the coefficients of f to those of f(Z*x)
  Zv = [2 3; 1 2]
  M = representation_matrix(representation_on_forms(G, 3))
  MZ = map_entries(e -> evaluate(e, QQ.(vec(Zv))), M)
  Rt, t = polynomial_ring(QQ, :t => 1:2)
  c = [1, 2, 3, 4]
  f = sum(c[i + 1] * binomial(3, i) * t[1]^(3 - i) * t[2]^i for i in 0:3)
  g = evaluate(f, Zv * t)
  @test MZ * c == [coeff(g, [3 - i, i])//binomial(3, i) for i in 0:3]

  #forms of a degree whose factorial exceeds the range of Int
  M = representation_matrix(representation_on_forms(G, 21))
  @test M[1, 2] == 21 * z[1, 1]^20 * z[2, 1]
  @test Oscar.multinomial_coefficient(30, [10, 10, 10]) == 5550996791340

  #direct sum of representation_linearly_reductive_group
  D = Oscar._direct_sum(rep1, rep2)
  @test ncols(representation_matrix(D)) == 6
  T1 = Oscar._tensor_product([rep1, rep2])
  T2 = Oscar._tensor_product(rep1, rep2)
  @test representation_matrix(T1) == representation_matrix(T2)
  @test ncols(representation_matrix(T1)) == 9

  #ternary cubics
  T, X = graded_polynomial_ring(QQ, :X => 1:10)
  g = linearly_reductive_group(:SL, 3, QQ)
  rep3 = representation_on_forms(g, 3)
  R_rep3 = invariant_ring(T, rep3)
  f =
    X[1] * X[4] * X[8] * X[10] - X[1] * X[4] * X[9]^2 - X[1] * X[5] * X[7] * X[10] +
    X[1] * X[5] * X[8] * X[9] + X[1] * X[6] * X[7] * X[9] - X[1] * X[6] * X[8]^2 -
    X[2]^2 * X[8] * X[10] + X[2]^2 * X[9]^2 + X[2] * X[3] * X[7] * X[10] -
    X[2] * X[3] * X[8] * X[9] + X[2] * X[4] * X[5] * X[10] - X[2] * X[4] * X[6] * X[9] -
    2 * X[2] * X[5]^2 * X[9] + 3 * X[2] * X[5] * X[6] * X[8] - X[2] * X[6]^2 * X[7] -
    X[3]^2 * X[7] * X[9] + X[3]^2 * X[8]^2 - X[3] * X[4]^2 * X[10] +
    3 * X[3] * X[4] * X[5] * X[9] - X[3] * X[4] * X[6] * X[8] - 2 * X[3] * X[5]^2 * X[8] +
    X[3] * X[5] * X[6] * X[7] + X[4]^2 * X[6]^2 - 2 * X[4] * X[5]^2 * X[6] + X[5]^4
  @test reynolds_operator(R_rep3, f) == f

  #tori
  #example in the derksen book
  T = torus_group(QQ, 1)
  r = representation_from_weights(T, [-3, -1, 1, 2])
  @test sprint(show, r) == "Representation of degree 4 of torus group"
  @test degree(r) == vector_space_dim(r) == 4
  @test_throws ArgumentError representation_from_weights(T, [-3 -1 1 2])
  @test weights(representation_from_weights(T, ZZ.([-3, -1, 1, 2]))) == weights(r)
  @test_throws ArgumentError representation_from_weights(torus_group(QQ, 2), [1, -1])
  I = invariant_ring(r)
  f = fundamental_invariants(I)
  @test length(f) == 6

  #several variables of the same weight
  I = invariant_ring(representation_from_weights(T, [1, 1, -1, -1]))
  X = gens(polynomial_ring(I))
  @test issetequal(fundamental_invariants(I), [X[i] * X[j] for i in 1:2 for j in 3:4])
  I = invariant_ring(representation_from_weights(T, [1, -1, 1, -1]))
  X = gens(polynomial_ring(I))
  @test issetequal(fundamental_invariants(I), [X[i] * X[j] for i in (1, 3) for j in (2, 4)])
  I = invariant_ring(representation_from_weights(T, [1, 1, -1, -1, -1]))
  X = gens(polynomial_ring(I))
  @test issetequal(fundamental_invariants(I), [X[i] * X[j] for i in 1:2 for j in 3:5])
  I = invariant_ring(representation_from_weights(T, [2, 2, -1, -1]))
  X = gens(polynomial_ring(I))
  @test issetequal(
    fundamental_invariants(I), [X[i] * X[j] * X[k] for i in 1:2 for j in 3:4 for k in j:4]
  )
  I = invariant_ring(representation_from_weights(T, [0, 0]))
  @test issetequal(fundamental_invariants(I), gens(polynomial_ring(I)))

  #all weights of the same sign: only the constants are invariant
  I = invariant_ring(representation_from_weights(T, [1, 2]))
  @test is_empty(fundamental_invariants(I))
  X = gens(polynomial_ring(I))
  @test reynolds_operator(I, X[1] * X[2] + 7) == 7
  @test_throws ArgumentError reynolds_operator(I, one(polynomial_ring(QQ, 3)[1]))
  @test is_zero(null_cone_ideal(I))
  @test base_ring(null_cone_ideal(I)) === polynomial_ring(I)

  #example from Macaulay2
  T = torus_group(QQ, 2)
  r = representation_from_weights(T, [1 0; 0 1; -1 -1; -1 1])
  I = invariant_ring(r)
  R = polynomial_ring(I)
  X = gens(R)
  f = fundamental_invariants(I)
  @test f == [X[1] * X[2] * X[3], X[1]^2 * X[3] * X[4]]

  #the Reynolds operator keeps the invariant terms, Hilbert's ideal is generated by
  #the fundamental invariants, and the Hilbert series counts invariant monomials
  for (m, W) in [(1, [-3, -1, 1, 2]), (2, [1 0; 0 1; -1 -1; -1 1])]
    r = representation_from_weights(torus_group(QQ, m), W)
    RT = invariant_ring(r)
    X = gens(polynomial_ring(RT))
    #t = (2, 3) multiplies x^a by 2^b[1]*3^b[2], where b is the weight of x^a
    t = [prod(QQ(p)^w for (p, w) in zip([2, 3], v)) for v in weights(r)]
    act = h -> evaluate(h, t .* X)
    N = 8
    f = (1 + sum(X))^N
    g = reynolds_operator(RT, f)
    @test act(g) == g
    @test all(u -> act(u) != u, monomials(f - g))
    FT = fundamental_invariants(RT)
    @test all(h -> reynolds_operator(RT, h) == h, FT)
    @test null_cone_ideal(RT) == ideal(FT)
    @test base_ring(null_cone_ideal(RT)) === polynomial_ring(RT)
    p, q = hilbert_series(RT)
    s = gen(parent(p))
    c = sum(s^total_degree(u) for u in monomials(g))
    @test is_zero(truncate(c * evaluate(q) - p, N + 1))
  end

  #another example, with affine algebra computation
  T = torus_group(QQ, 2)
  r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1])
  RT = invariant_ring(r)
  @test coefficient_ring(RT) === QQ
  A, _ = affine_algebra(RT)
  @test ngens(modulus(A)) == 1

  #no reynolds operator required/available:

  #SL(2) over symmetric forms of degree 2
  g = linearly_reductive_group(:SL, 2, QQ)
  r = representation_on_forms(g, 2)
  M = representation_matrix(r)
  ringg = parent(M[1, 1])
  z = gens(ringg)
  f = z[1] * z[4] - z[2] * z[3] - 1
  G = linearly_reductive_group(ideal([f]))
  @test_throws ArgumentError representation_on_forms(G, 2)
  @test_throws ArgumentError representation_reductive_group(G)
  @test_throws ArgumentError Oscar.natural_representation(G)
  R = representation_reductive_group(G, M)
  RG = invariant_ring(R)
  F = fundamental_invariants(RG)
  @test length(F) == 1
  X = polynomial_ring(RG)
  g = gens(X)
  @test F[1] == g[1] * g[3] - g[2]^2

  #a torus acting with weights 1, 1 has no invariants of positive degree
  ringg, z = polynomial_ring(QQ, :z => 1:2)
  G = linearly_reductive_group(ideal([z[1] * z[2] - 1]))
  @test base_ring(G) === QQ
  Zz, w = polynomial_ring(ZZ, :w => 1:2)
  @test_throws ArgumentError linearly_reductive_group(ideal([w[1] * w[2] - 1]))
  R = representation_reductive_group(G, diagonal_matrix(ringg, [z[1], z[1]]))
  @test is_empty(fundamental_invariants(invariant_ring(R)))

  #tori given by their ideal agree with the torus algorithm
  for w in [[-3, -1, 1, 2], [1, 1, -1, -1], [2, -1, -1]]
    IT = invariant_ring(representation_from_weights(torus_group(QQ, 1), w))
    M = diagonal_matrix(ringg, [a >= 0 ? z[1]^a : z[2]^-a for a in w])
    RG = invariant_ring(polynomial_ring(IT), representation_reductive_group(G, M))
    F = fundamental_invariants(RG)
    FT = fundamental_invariants(IT)
    @test length(F) == length(FT)
    @test all(f -> subalgebra_membership_homogeneous(f, FT)[1], F)
  end
  ringg, z = polynomial_ring(QQ, :z => 1:4)
  G = linearly_reductive_group(ideal([z[1] * z[2] - 1, z[3] * z[4] - 1]))
  M = diagonal_matrix(ringg, [z[1], z[3], z[2] * z[4], z[2] * z[3]])
  IT = invariant_ring(
    representation_from_weights(torus_group(QQ, 2), [1 0; 0 1; -1 -1; -1 1])
  )
  RG = invariant_ring(polynomial_ring(IT), representation_reductive_group(G, M))
  @test_throws ArgumentError fundamental_invariants(RG, :reynolds)
  @test issetequal(fundamental_invariants(RG), fundamental_invariants(IT))

  #SLm with and without the Reynolds operator
  g = linearly_reductive_group(:SL, 2, QQ)
  r4 = representation_on_forms(g, 4)
  for r in [
    representation_on_forms(g, 3),
    r4,
    Oscar._direct_sum(representation_on_forms(g, 2), representation_on_forms(g, 2)),
    representation_on_forms(linearly_reductive_group(:SL, 3, QQ), 2),
  ]
    rg = invariant_ring(r)
    FF = fundamental_invariants(rg, :reynolds)
    @test all(h -> is_invariant(h, r), FF)
    F = fundamental_invariants(invariant_ring(polynomial_ring(rg), r), :linear_algebra)
    @test length(F) == length(FF)
    @test all(f -> subalgebra_membership_homogeneous(f, FF)[1], F)
  end
  @test length(fundamental_invariants(invariant_ring(r4))) == 2
  @test_throws ArgumentError fundamental_invariants(invariant_ring(r4), :king)

  #binary quartics: the invariants of degrees 2 and 3 are algebraically independent
  p, q = hilbert_series(invariant_ring(r4))
  t = gen(parent(p))
  @test p * (1 - t^2) * (1 - t^3) == evaluate(q)
end

@testset "Supertype of invariant rings" begin
  #the invariant rings of finite groups, linearly reductive groups and tori share a supertype
  @test Oscar.FinGroupInvarRing <: Oscar.InvarRing
  @test Oscar.RedGroupInvarRing <: Oscar.InvarRing
  @test Oscar.TorGroupInvarRing <: Oscar.InvarRing
end
