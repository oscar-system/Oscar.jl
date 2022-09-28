@testset "fraction fields of varieties" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  @test is_irreducible(Spec(R))
  @test is_irreducible(Spec(R, ideal(R, x)))
  @test !is_irreducible(Spec(R, ideal(R, x*y)))
  @test is_irreducible(Spec(Localization(R, units_of(R))[1]))
  @test !is_irreducible(Spec(R, ideal(R, x*y), units_of(R)))

  P = projective_space(QQ, 2)
  S = homogeneous_poly_ring(P)
  C = subscheme(P, ideal(S, S[1]*S[2]-S[3]^2))
  Ccov = as_covered_scheme(C)

  KK = VarietyFunctionField(Ccov)
  U2 = patches(Ccov)[2]
  a = 5*gens(OO(U2))[1]*gens(OO(U2))[2]

  b = KK(a, one(a))
  @test b^2 - b == KK(a^2-a, one(a))
end

@testset "fraction fields of varieties II" begin
  # We construct by hand the projective bundle
  # P(O_{P^1}(4)+O_{P^1}(6)+O_{P^1}(1))
  kk = GF(29)

  # Set up the base ℙ¹ with coordinates s and t
  R, (s,t) = PolynomialRing(kk, ["s", "t"])
  S, _ = grade(R, [1, 1])

  base_P1 = ProjectiveScheme(S)

  # split this into the standard covering
  base_covering = standard_covering(base_P1)

  A1s = patches(base_covering)[1]
  A1t = patches(base_covering)[2]

  # Set up relative projective space of relative dimension 2 
  # over both base patches
  P2_s = projective_space(OO(A1s), ["xs", "ys", "zs"])

  Cs = standard_covering(P2_s)

  P2_t = projective_space(OO(A1t), ["xt", "yt", "zt"])

  Ct = standard_covering(P2_t)

  # Join the resulting schemes in a disjoint union with two 
  # components
  C = disjoint_union(Cs, Ct)

  # Manually glue two (dense) patches of the two components
  X = Cs[3]
  Y = Ct[3]
  x = gens(base_ring(OO(X)))
  y = gens(base_ring(OO(Y)))
  f = maximal_extension(X, Y, [x[1]//(x[3])^4, x[2]//(x[3])^6, 1//x[3]])
  g = maximal_extension(Y, X, [y[1]//(y[3])^4, y[2]//(y[3])^6, 1//y[3]])
  add_glueing!(C, Glueing(X, Y, restriction(f, domain(f), domain(g)), restriction(g, domain(g), domain(f))))

  # Extend the glueing to the whole covered scheme
  fill_transitions!(C)

  X = CoveredScheme(C)

  KK = VarietyFunctionField(X)

  U = representative_patch(KK)
  V = C[5]
  R = ambient_ring(V)
  x = gens(R)
  f = x[1]^2 - 2*x[2]^5*x[3]^3
  g = 4*x[3]^2 - 5*x[2]
  @test KK(f, g) + KK(f, g) == KK(2*f//g)
  h = KK(f, g)
  @test h[V] == f//g
  K = FractionField(R)
  @test K(h) == f//g

  @test KK(f+2*g-5, g) + KK(f+2*g-5, g) == KK(2*(f+2*g-5)//g)
  h = KK(f+2*g-5, g)
  @test h[V] == (f+2*g-5)//g
  K = FractionField(R)
  @test K(h) == (f+2*g-5)//g
end
