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
  add_glueing!(C, Glueing(X, Y, restrict(f, domain(f), domain(g)), restrict(g, domain(g), domain(f))))

  # Extend the glueing to the whole covered scheme
  fill_transitions!(C)

  X = CoveredScheme(C)

  U = C[1]
  x = gens(ambient_coordinate_ring(U))
# I = IdealSheaf(X, U, OO(U).([x[1]-1, x[2]-2, x[3]-3]))
# J = IdealSheaf(X, U, OO(U).([x[1]-5, x[2]-1, x[3]]))
  I = IdealSheaf(X, U, OO(U).([x[1]-1]))
  J = IdealSheaf(X, U, OO(U).([x[2]-5]))
  D = WeilDivisor(I)
  E = WeilDivisor(J)
  @test D + 2*E == D + E + E

  KK = VarietyFunctionField(X)
  U = representative_patch(KK)
  V = C[3]
  R = ambient_coordinate_ring(V)
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

  # create an elliptically fibered K3
  Ut = Ct[3]
  x, y, t = coordinates(Ut)
  ft = y^2 - (x^3 + 21*x + (28*t^7+18))
  I = IdealSheaf(X, Ut, [ft])
  adeK3 = subscheme(I)
  @test dim(singular_locus(adeK3)[1]) == 0

  x,y,t = coordinates(adeK3[1][6])
  # ideal defining a section of the fibration
  P = [(5*t^8 + 20*t^7 + 2*t^6 + 23*t^5 + 20*t^3 + 11*t^2 + 3*t + 13) - x*(t^4 + 9*t^3 + 5*t^2 + 22*t + 16),
  (26*t^12 + 11*t^11 + 15*t^10 + 8*t^9 + 20*t^8 + 25*t^7 + 16*t^6 + 3*t^5 + 10*t^4 + 15*t^3 + 4*t^2 + 28*t + 28) - y*(t^6 + 28*t^5 + 27*t^4 + 23*t^3 + 8*t^2 + 13*t + 23)]
  # the following computation dies ... but it should not.
  D = IdealSheaf(adeK3, adeK3[1][6], P)
  Dscheme = subscheme(D)

  # an interesting rational function
  K = function_field(adeK3)
  x,y,t= ambient_coordinates(adeK3[1][6])
  phi = K(-8*t^8 + 4*t^7 + 6*t^5 + 4*x*t^3 + 3*t^4 - 13*x*t^2 - 7*y*t^2 - 12*t^3 - 12*x*t + 12*y*t - 14*t^2 - 14*x - y - 10*t + 10)//K(6*t^8 - 5*t^7 + 14*t^6 - 7*x*t^4 - 13*t^5 - 5*x*t^3 - 6*x*t^2 - 5*t^3 - 9*x*t - 10*t^2 + 4*x - 8*t + 4)
  @test Oscar.order_on_divisor(phi, D) == -1
end

@testset "orders on divisors" begin
  kk = QQ
  R, (s,t) = PolynomialRing(kk, ["s", "t"])
  X = Spec(R)
  Xc = CoveredScheme(X)
  KK = VarietyFunctionField(Xc)
  f = s^2 + t^2-1
  I = IdealSheaf(Xc, X, [f])
  F = KK(f^70)
  @test order_on_divisor(F, I) == 70
end

@testset "linear systems" begin
  P2 = projective_space(QQ, 2)
  S = ambient_coordinate_ring(P2)
  X = covered_scheme(P2)
  I = IdealSheaf(P2, [S[1]])
  D = WeilDivisor(I)

  KK = function_field(X)
  R = ambient_coordinate_ring(representative_patch(KK))
  x = gens(R)
  @test in_linear_system(KK(x[1]), D)
  @test !in_linear_system(KK(x[1]^2), D)
  @test in_linear_system(KK(x[1]^2), 2*D)
  # Not running at the moment; work in progress
  #@test !in_linear_system(KK(x[1], x[2]), D)

  L = LinearSystem(KK.([1, x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2]), 2*D)
  H = S[1]+S[2]+S[3]
  P = IdealSheaf(P2, [H])
  @test ngens(subsystem(L, P, 1)[1]) == 3
  @test ngens(subsystem(L, P, 2)[1]) == 1
  @test ngens(subsystem(L, P, 3)[1]) == 0
end

@testset "WeilDivisor" begin
  P2 = projective_space(QQ, 2)
  S = ambient_coordinate_ring(P2)
  (x, y, z) = gens(S)
  I = ideal(S, x^3 + y^3 + z^3)
  C = subscheme(P2, I)
  X = covered_scheme(C)
  H = ideal(S, [x, y+z])
  J = ideal_sheaf(C, H)

  coeff_dict = IdDict{IdealSheaf, fmpz}()
  coeff_dict[J] = ZZ(3)
  D = WeilDivisor(X, ZZ, coeff_dict)
  @test D == 3*weil_divisor(J)
  @test sprint(show, D-D) isa String
  @test (ZZ(2)*D == D + J + J + J)
  @test !(D-D == 4*D)
  KK = function_field(X)
  U = X[1][2]
  u, v = gens(OO(U))
  f = KK(u, v) 
  @test !in_linear_system(f, D)
end
