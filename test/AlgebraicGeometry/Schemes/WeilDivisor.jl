@testset "Ideal sheaves and Weil divisors" begin
  kk = GF(29)

  # Set up the base ℙ¹ with coordinates s and t
  S, (s, t) = graded_polynomial_ring(kk, [:s, :t])

  base_P1 = proj(S)

  # split this into the standard covering
  bc = standard_covering(base_P1)

  A1s = patches(bc)[1]
  A1t = patches(bc)[2]

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
  x = gens(OO(X))
  y = gens(OO(Y))
  f = maximal_extension(X, Y, [x[1]//(x[3])^4, x[2]//(x[3])^6, 1//x[3]])
  g = maximal_extension(Y, X, [y[1]//(y[3])^4, y[2]//(y[3])^6, 1//y[3]])
  add_gluing!(C, Gluing(X, Y, restrict(f, domain(f), domain(g)), restrict(g, domain(g), domain(f))))

  # Extend the gluing to the whole covered scheme
  Oscar.fill_transitions!(C)

  X = CoveredScheme(C)

  U = C[1]
  x = gens(ambient_coordinate_ring(U))
# I = IdealSheaf(X, U, OO(U).([x[1]-1, x[2]-2, x[3]-3]))
# J = IdealSheaf(X, U, OO(U).([x[1]-5, x[2]-1, x[3]]))
  I = IdealSheaf(X, U, OO(U).([x[1]-1]))
  J = IdealSheaf(X, U, OO(U).([x[2]-5]))
  Oscar.maximal_associated_points(I)
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
  K = fraction_field(R)
  @test K(h) == f//g

  @test KK(f+2*g-5, g) + KK(f+2*g-5, g) == KK(2*(f+2*g-5)//g)
  h = KK(f+2*g-5, g)
  @test h[V] == (f+2*g-5)//g
  K = fraction_field(R)
  @test K(h) == (f+2*g-5)//g

  # create an elliptically fibered K3
  Ut = Ct[3]
  x, y, t = coordinates(Ut)
  ft = y^2 - (x^3 + 21*x + (28*t^7+18))
  I = IdealSheaf(X, Ut, [ft])
  adeK3, inc_adeK3 = sub(I)
  @test dim(singular_locus(adeK3)[1]) == 0

  weier_chart = first([U for U in affine_charts(adeK3) if codomain(covering_morphism(inc_adeK3)[U]) === X[1][6]]) # Order of charts is random due to use of dictionaries in the constructor
  x,y,t = coordinates(weier_chart)
  # ideal defining a section of the fibration
  P = [(5*t^8 + 20*t^7 + 2*t^6 + 23*t^5 + 20*t^3 + 11*t^2 + 3*t + 13) - x*(t^4 + 9*t^3 + 5*t^2 + 22*t + 16),
  (26*t^12 + 11*t^11 + 15*t^10 + 8*t^9 + 20*t^8 + 25*t^7 + 16*t^6 + 3*t^5 + 10*t^4 + 15*t^3 + 4*t^2 + 28*t + 28) - y*(t^6 + 28*t^5 + 27*t^4 + 23*t^3 + 8*t^2 + 13*t + 23)]
  # the following computation dies ... but it should not.
  D = IdealSheaf(adeK3, weier_chart, P)
  Dscheme, inc_Dscheme = sub(D)

  # an interesting rational function
  K = function_field(adeK3)
  x,y,t= ambient_coordinates(weier_chart)
  phi = K(-8*t^8 + 4*t^7 + 6*t^5 + 4*x*t^3 + 3*t^4 - 13*x*t^2 - 7*y*t^2 - 12*t^3 - 12*x*t + 12*y*t - 14*t^2 - 14*x - y - 10*t + 10)//K(6*t^8 - 5*t^7 + 14*t^6 - 7*x*t^4 - 13*t^5 - 5*x*t^3 - 6*x*t^2 - 5*t^3 - 9*x*t - 10*t^2 + 4*x - 8*t + 4)
  @test Oscar.order_of_vanishing(phi, D, check=false) == -1

  other_chart = first([U for U in affine_charts(adeK3) if codomain(covering_morphism(inc_adeK3)[U]) === X[1][2]]) # Order of charts is random due to use of dictionaries in the constructor
  (x,z,t) = coordinates(other_chart)
  o = weil_divisor(ideal_sheaf(adeK3, other_chart, [z,x]), check=false)
  @test order_of_vanishing(K(z), o, check=false) == 3
  @test order_of_vanishing(K(x), o, check=false) == 1
end

@testset "orders on divisors" begin
  kk = QQ
  R, (s,t) = polynomial_ring(kk, [:s, :t])
  X = spec(R)
  Xc = CoveredScheme(X)
  KK = VarietyFunctionField(Xc)
  f = s^2 + t^2-1
  I = IdealSheaf(Xc, X, [f])
  F = KK(f^70)
  @test order_of_vanishing(F, I) == 70
end

@testset "linear systems" begin
  P2 = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(P2)
  X = covered_scheme(P2)
  I = IdealSheaf(P2, [S[1]])
  D = WeilDivisor(I)

  KK = function_field(X)
  R = ambient_coordinate_ring(representative_patch(KK))
  x = gens(R)
  @test is_in_linear_system(KK(x[1]), D)
  @test !is_in_linear_system(KK(x[1]^2), D)
  @test is_in_linear_system(KK(x[1]^2), 2*D)
  @test !is_in_linear_system(KK(x[1], x[2]), D)

  L = LinearSystem(KK.([1, x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2]), 2*D)
  H = S[1]+S[2]+S[3]
  P = IdealSheaf(P2, [H])
  @test ngens(Oscar._subsystem(L, P, 1)[1]) == 3
  @test ngens(Oscar._subsystem(L, P, 2)[1]) == 1
  @test ngens(Oscar._subsystem(L, P, 3)[1]) == 0
end

@testset "WeilDivisor" begin
  P2 = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(P2)
  (x, y, z) = gens(S)
  I = ideal(S, x^3 + y^3 + z^3)
  C = subscheme(P2, I)
  SC = homogeneous_coordinate_ring(C)
  (x, y, z) = gens(SC)
  X = covered_scheme(C)
  H = ideal(SC, SC.([x, y+z]))
  J = ideal_sheaf(C, H)

  coeff_dict = IdDict{IdealSheaf, ZZRingElem}()
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
  @test !is_in_linear_system(f, D)
end

@testset "intersections of weil divisors on surfaces" begin
  P = projective_space(QQ, 2)
  X = covered_scheme(P)
  S = homogeneous_coordinate_ring(P)
  x = gens(S)
  I1 = IdealSheaf(P, x[1]^3 + x[2]^3 + x[3]^3)
  I2 = IdealSheaf(P, x[2])
  I3 = IdealSheaf(P, x[1] + x[2] + x[3])

  D1 = 7 * weil_divisor(I1)
  D2 = weil_divisor(I2)
  D3 = weil_divisor(I3)
  @test intersect(D1, D2) == intersect(D1, D3) == 21
  @test intersect(D2, D3) == 1
end

@testset "decomposition" begin
  P3 = projective_space(QQ, 3)
  S = homogeneous_coordinate_ring(P3)
  (x, y, z, w) = gens(S)
  I = ideal(S, x^2*y^3)
  II = ideal_sheaf(P3, I)
  X = covered_scheme(P3)
  D = weil_divisor(II)
  E = Oscar.irreducible_decomposition(D)
  @test length(keys(Oscar.coefficient_dict(E))) == 2
  @test 2*one(coefficient_ring(E)) in values(Oscar.coefficient_dict(E))
  @test 3*one(coefficient_ring(E)) in values(Oscar.coefficient_dict(E))
end

@testset "intersection numbers on surfaces" begin
  P3 = projective_space(QQ, 3)
  S = homogeneous_coordinate_ring(P3)
  (x,y, z, w) = gens(S)
  I = ideal(S, [x^4+y^4+z^4+w^4])
  II = ideal_sheaf(P3, I)
  P = covered_scheme(P3)
  inc = Oscar.CoveredClosedEmbedding(covered_scheme(P3), II)
  X = domain(inc)
  C1 = EffectiveCartierDivisor(ideal_sheaf(P3, [x+y+z+w]))
  C2 = EffectiveCartierDivisor(ideal_sheaf(P3, [x^2*y + y^2*z + z^2*w + w^2*x]))
  C1 = pullback(inc)(C1)
  C2 = pullback(inc)(C2)
  d = intersect(weil_divisor(C1), weil_divisor(C2))
  pts = Oscar.irreducible_decomposition(intersect(C1, C2))
  @test integral(pts) == d
end
