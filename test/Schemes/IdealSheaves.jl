@testset "Ideal sheaves and Weil divisors" begin
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

  U = C[1]
  x = gens(ambient_ring(U))
  I = IdealSheaf(X, [x[1]-1, x[2]-2, x[3]-3])
  set_name!(I, "I")
  J = IdealSheaf(X, [x[1]-5, x[2]-1, x[3]])
  set_name!(J, "J")
  D = WeilDivisor(I)
  E = WeilDivisor(J)
  @test D + 2*E == D + E + E

  KK = VarietyFunctionField(X)
  U = representative_patch(KK)
  V = C[3]
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
  #order_on_divisor(one(KK), I)
end
