@testset "Ideal sheaves" begin
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
  add_glueing!(C, Glueing(X, Y, restrict(f, domain(f), domain(g)), restrict(g, domain(g), domain(f))))

  # Extend the glueing to the whole covered scheme
  fill_transitions!(C)

  X = CoveredScheme(C)

  U = C[1]
  x = gens(ambient_ring(U))
  I = IdealSheaf(X, U, OO(U).([x[1]-1, x[2]-2, x[3]-3]))
  J = IdealSheaf(X, U, OO(U).([x[1]-5, x[2]-1, x[3]]))

  @test I+J == IdealSheaf(X, U, [one(OO(U))])
  K = I*J
end

@testset "ideal sheaves II" begin
  IP2 = projective_space(QQ, 2)
  X = covered_scheme(IP2)
  S = ambient_ring(IP2)
  (u,v,w) = gens(S)
  Ihom = ideal(S, u^2 - v*w)
  I = IdealSheaf(IP2, Ihom)
  @test is_prime(I)
  I2 = IdealSheaf(IP2, gens(Ihom))
  I3 = IdealSheaf(IP2, gens(Ihom)[1])
  @test I == I2 == I3
  U = patches(default_covering(X))
  @test I(U[1]) isa Ideal 
  @test I(U[2]) isa Ideal 
  @test I(U[3]) isa Ideal 
  V = PrincipalOpenSubset(U[1], gens(OO(U[1]))[1])
  rho = I(U[1], V)
  @test I(V) == ideal(OO(V), rho.(gens(I(U[1]))))
  Y = subscheme(I)
  simplify!(I)

  # run the check block in extend!
  ID = IdDict{AbsSpec, Ideal}()
  ID[X[1][1]] = I(X[1][1])
  J = IdealSheaf(X, extend!(default_covering(X), ID), check=true)
  ID[X[1][1]] = I(X[1][1])
  J2 = IdealSheaf(X, ID)
  @test J == J2
  @test sprint(show, J) isa String

  @test issubset(IdealSheaf(X), I) # Whether the zero ideal sheaf is a subset of I
end
