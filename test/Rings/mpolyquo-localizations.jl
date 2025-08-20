@testset "mpolyquo-localizations.jl" begin
  R, v = QQ[:x, :y, :u, :v]
  x = v[1]
  y = v[2] 
  u = v[3]
  v = v[4]
  f = x*v-y*u
  I = ideal(R, f)
  Q, p = quo(R, I)
  S = Oscar.MPolyComplementOfKPointIdeal(R, [QQ(1), QQ(0), QQ(1), QQ(0)])
  L, _ = localization(Q, S)
  a = L(x)
  b = L(y)
  c = L(u)
  d = L(v)
  
  kk= QQ
  R, v = kk[:x, :y]
  x = v[1]
  y = v[2] 
  f = (x^2 + y^2)
  T = Oscar.MPolyComplementOfKPointIdeal(R, [kk(0), kk(0)])
  I = ideal(R, f)
  V = Oscar.MPolyQuoLocRing(R, I, T)

  S = R
  U = Oscar.MPolyPowersOfElement(S, [f-1])
  J = ideal(S, zero(S))
  W = Oscar.MPolyQuoLocRing(S, J, U)

  h = Oscar.MPolyQuoLocalizedRingHom(W, V, [x//(y-1), y//(x-5)])
  J1 = ideal(V, [x*(x-1), y*(y-3)])
  J2 = ideal(W, [x, y])
  @test preimage(h, J1) == J2
  
  ### second round of tests
  #kk = GF(101)
  ⊂ = issubset
  kk = QQ
  R, (x,y) = kk[:x, :y]

  f = x^2 + y^2-2
  S = Oscar.MPolyPowersOfElement(x-1)
  T = Oscar.MPolyComplementOfKPointIdeal(R, [1,1])
  V = Oscar.MPolyComplementOfPrimeIdeal(ideal(R, f))
  ⊂ = issubset
  @test S ⊂ V
  @test !(V ⊂ S)
  @test !(T ⊂ V)
  @test (V ⊂ T)
  @test !(Oscar.MPolyComplementOfPrimeIdeal(ideal(R, f-1)) ⊂ T)
  @test S ⊂ Oscar.MPolyComplementOfPrimeIdeal(ideal(R, f-1))
  @test !(Oscar.MPolyPowersOfElement(f) ⊂ V)
  @test Oscar.MPolyPowersOfElement(x-1) ⊂ Oscar.MPolyComplementOfKPointIdeal(R, [0,0])
  @test Oscar.MPolyPowersOfElement(x-1) * Oscar.MPolyComplementOfKPointIdeal(R, [0,0]) ⊂ Oscar.MPolyComplementOfKPointIdeal(R, [0,0])
  @test !(Oscar.MPolyPowersOfElement(x) ⊂ Oscar.MPolyComplementOfKPointIdeal(R, [0,0]))
  @test T*T == T

  U = S*T
  @test U[1] ⊂ S || U[1] ⊂ T
  @test U[2] ⊂ S || U[2] ⊂ T
  @test S ⊂ U 
  @test T ⊂ U
  @test S*U == U
  @test T*U == U 
  @test U*U == U
  g = rand(S, 0:3, 1:5, 2:8)
  g = rand(T, 0:3, 1:5, 2:8)
  g = rand(U, 0:3, 1:5, 2:8)
  W, _ = localization(U)
  localization(W, S)
  @test base_ring(W) == R
  @test inverted_set(W) == U
  L, _ = quo(W, ideal(W, f))
  @test issubset(modulus(underlying_quotient(L)), saturated_ideal(modulus(L)))
  @test gens(L) == L.(gens(R))
  @test x//(y*(x-1)) in L

  I = ideal(L, (x-1)*(y-1))
  @test one(L) in I
  @test is_unit(L(y-1))

  h = x^4+23*x*y^3-15
  Q, _ = quo(R, f)
  T = Oscar.MPolyPowersOfElement(h^3)
  W, _ = localization(Q, T)
  @test x//(h+3*f) in W
  @test W(x//(h+3*f)) == W(x//h)
  g = W.([rand(R, 0:5, 0:2, 0:1) for i in 1:10])
  @test isone(W(h)*inv(W(h)))
  g = [W(8*x*h+4), W(x, 45*h^2)]
  @test one(localized_ring(W)) in ideal(W, g)

  @test one(W) == dot(Oscar.write_as_linear_combination(one(W), g), g)


  h = (x+5)*(x^2+10*y)+(y-7)*(y^2-3*x)
  Q, _ = quo(R, h)
  T = Oscar.MPolyComplementOfKPointIdeal(R, [-5, 7])
  W, _ = localization(Q, T)
  @test x//(y) in W
  @test x//(y+h) in W
  g = [W(h + (x+5) - 9, y+24*x^3-8)]
  (d, a) = Oscar.bring_to_common_denominator([one(W), g[1]])
  @test W(a[1], d) == one(W)
  @test W(a[2]*lifted_numerator(g[1]), d) == g[1]
  @test (one(localized_ring(W)) in ideal(W, g))
  @test one(W) == dot(Oscar.write_as_linear_combination(one(W), g), g) 
end

@testset "prime ideals in quotient rings" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, [x^2-y^2])
  U = powers_of_element(y)
  A, = quo(R, I)
  L, = localization(A, U)
  J = ideal(R, [y*(x-1), y*(y+1)])
  @test !is_prime(J)
  @test is_prime(L(J))
end

@testset "additional hom constructors" begin
  R, (x, y) = ZZ[:x, :y]
  I = ideal(R, [x^2-y^2])
  U = powers_of_element(y)
  A, = quo(R, I)
  L, = localization(A, U)
  phi = hom(L, L, gens(L))
end
  
@testset "fractions and divexact: Issue #1885" begin
  R, (x, y) = ZZ[:x, :y]
  I = ideal(R, [x^2-y^2])
  U = powers_of_element(y)
  A, = quo(R, I)
  L, = localization(A, U)

  @test L(x)/L(y) == L(x//y)
  @test L(y)/L(x) == L(y//x)
  @test_throws ErrorException L(y)/L(x-1)
end

@testset "associated primes (quo and localized)" begin
  # define the rings
  R,(x,y,z,w) = QQ[:x, :y, :z, :w]
  Q1 = ideal(R,[x^2+y^2-1])
  Q2 = ideal(R,[x*y-z*w])
  RQ1,phiQ1 = quo(R,Q1)
  RQ2,phiQ2 = quo(R,Q2)
  T1 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0,0])
  f = x+y+z+w-1
  T2 = Oscar.MPolyPowersOfElement(f)
  RL1,phiL1 = localization(R,T1)
  RL2,phiL2 = localization(R,T2)
  RQ1L1, phiQ1L1 = localization(RQ1,T1)
  RQ1L2, phiQ1L2 = localization(RQ1,T2)
  RQ2L1, phiQ2L1 = localization(RQ2,T1)
  RQ2L2, phiQ2L2 = localization(RQ2,T2)

  # tests for MPolyQuoRing
  I1 = ideal(R,[x*z*(w-1),y*z*(w-1)])
  I2 = ideal(R,[y^2*x*(x+y+w-1),z])
  LM = minimal_primes(phiQ1(I1))
  LP = primary_decomposition(phiQ1(I1))  ## coincides with min. ass. primes here
  @test length(LM) == 2
  @test length(LP) == 2
  @test intersect(LM) == phiQ1(radical(I1+Q1))
  @test intersect([a for (a,_) in LP]) == intersect(LM)
  LM = minimal_primes(phiQ2(I1))
  LP = primary_decomposition(phiQ2(I1)) 
  @test length(LM) == 4
  @test intersect(LM) == phiQ2(radical(I1+Q2))
  @test length(LP) == 5
  @test intersect([a for (a,_) in LP]) == phiQ2(I1+Q2)

  # tests for MPolyLocRing
  LM = minimal_primes(phiL1(I1))   
  LP = primary_decomposition(phiL1(I1)) ## coincides with min. ass. primes here
  @test length(LM) == 2
  @test length(LP) == 2
  @test intersect(LM) == phiL1(I1)
  @test intersect([a for (a,_) in LP]) == phiL1(I1)
  LM = minimal_primes(phiL2(I2))
  LP = primary_decomposition(phiL2(I2))
  @test intersect(LM) == ideal(RL2,RL2.([z,x*y]))
  @test intersect([a for (a,_) in LP]) == ideal(RL2,RL2.([z,x*y^2]))

  # tests for MPolyLocQuoRing
  LM = minimal_primes(phiQ1L1(phiQ1(I1))) 
  @test length(LM)==0
  LM = minimal_primes(phiQ2L1(phiQ2(I1)))
  LP = primary_decomposition(phiQ2L1(phiQ2(I1)))
  @test length(LM) == 3
  @test length(LP) == 4
  @test intersect(LM) == phiQ2L1(phiQ2(radical(I1+Q2)))
  @test intersect([a for (a,_) in LP]) == phiQ2L1(phiQ2(I1+Q2))
  @test length(minimal_primes(phiQ2L1(phiQ2(I2)))) == 2
  @test length(minimal_primes(phiQ2L2(phiQ2(I2)))) == 2
end

@testset "saturation (quo and localization)" begin
  # define the rings
  R,(x,y,z) = QQ[:x, :y, :z]
  Q1 = ideal(R,[x])
  Q2 = ideal(R,[z,x^2-y^2])
  RQ1,phiQ1 = quo(R,Q1)
  RQ2,phiQ2 = quo(R,Q2)
  T1 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  f = x-y
  T2 = Oscar.MPolyPowersOfElement(f)
  RL1,phiL1 = localization(R,T1)
  RL2,phiL2 = localization(R,T2)
  RQ1L1, phiQ1L1 = localization(RQ1,T1)
  RQ1L2, phiQ1L2 = localization(RQ1,T2)
  RQ2L1, phiQ2L1 = localization(RQ2,T1)
  RQ2L2, phiQ2L2 = localization(RQ2,T2)

  # the ideals
  I1 = ideal(R,[x^2*y+y^2*z+z^2*x])
  I2 = ideal(R,[x-y])
  I3 = ideal(R,[x,y,z])
  I4 = ideal(R,[x-1,y-1,z])
  I5 = ideal(R,[x,y+1,z-1])
  I6 = ideal(R,[y])

  # quotient ring tests
  I = phiQ1(I1)
  J = phiQ1(I6)
  @test saturation_with_index(I,J)[2] == 2
  K = phiQ1(I2*I3*I4*I5)
  @test saturation(K,J) == phiQ1(I5)
  I = phiQ2(I1)
  J = phiQ2(I6)
  @test saturation_with_index(I,J)[2] == 3
  K = phiQ2(I2*I3*I4*I5)
  @test saturation(K,J) == phiQ2(I2*I4)

  # localized ring tests
  I = phiL1(I1)
  J = phiL1(I6)
  @test saturation_with_index(I,J)[2] == 0
  K = phiL1(I2*I3*I4*I5)
  @test saturation(K,J) == phiL1(I2)
  I = phiL2(I1*I6^2)
  J = phiL2(I6)
  @test saturation_with_index(I,J)[2] == 2
  K = phiL2(I2*I3*I4*I5)
  @test saturation(K,J) == phiL2(I5)

  # localized quo ring tests
  I = phiQ1L1(phiQ1(I1))
  J = phiQ1L1(phiQ1(I6))
  @test saturation_with_index(I,J)[2] == 2
  K = phiQ1L1(phiQ1(I2*I3*I4*I5))
  @test saturation(K,J) == ideal(RQ1L1,one(RQ1L1))
  @test saturation_with_index(K,J)[2] == 2
  I = phiQ1L2(phiQ1(I1))
  J = phiQ1L2(phiQ1(I6))
  @test saturation(I,J) == phiQ1L2(phiQ1(ideal(R,[x,z])))
  K = phiQ1L2(phiQ1(I2*I3*I4*I5))
  @test saturation(K,J) == phiQ1L2(phiQ1(ideal(R,[z - 1, y + 1])))
  I = phiQ2L1(phiQ2(I1))
  J = phiQ2L1(phiQ2(I6))
  @test saturation_with_index(I,J)[2] == 3
  K = phiQ2L1(phiQ2(I2*I3*I4*I5))
  @test saturation(K,J) == phiQ2L1(phiQ2(ideal(R,[x-y])))
end

@testset "algebras as vector spaces" begin
  R, (x,y) = QQ[:x, :y]
  I = ideal(R, [x^3- 2, y^2+3*x])
  I = (x-8)^2*I
  L, _ = localization(R, powers_of_element(x-8))
  A, pr = quo(L, L(I))
  V, id = vector_space(QQ, A)
  @test vector_space_dim(V) == 6
  @test id.(gens(V)) == A.([x^2*y, x*y, y, x^2, x, one(x)])
  f = (x*3*y-4)^5
  f = A(f)
  @test id(preimage(id, f)) == f
  v = V[3] - 4*V[5]
  @test preimage(id, id(v)) == v

  R, (x,y) = QQ[:x, :y]
  I = ideal(R, [x^3, (y-1)^2+3*x])
  I = (x-8)^2*I
  L, _ = localization(R, complement_of_point_ideal(R, [0, 1]))
  A, pr = quo(L, L(I))
  V, id = vector_space(QQ, A)
  @test vector_space_dim(V) == 6
  f = (x*3*y-4)^5
  f = A(f)
  @test id(preimage(id, f)) == f
  v = V[3] - 4*V[5]
  @test preimage(id, id(v)) == v
end

@testset "minimal and small generating sets" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  IQ = ideal(R,[x-z])
  U1 = Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  U2 = Oscar.MPolyComplementOfKPointIdeal(R,[1,1,1])
  U3 = Oscar.MPolyPowersOfElement(y)
  Q1 = Oscar.MPolyQuoLocRing(R,IQ,U1)
  Q2 = Oscar.MPolyQuoLocRing(R,IQ,U2)
  Q3 = Oscar.MPolyQuoLocRing(R,IQ,U3)
  J1 = ideal(Q1,[x^2-y^2,y^2-z^2,x^2-z^2])
  @test length(minimal_generating_set(J1)) == 1
  @test length(small_generating_set(J1)) == 1
  @test ideal(Q1, small_generating_set(J1)) == ideal(Q1, minimal_generating_set(J1))
  @test ideal(Q1, minimal_generating_set(J1)) == J1

  J2 = ideal(Q2,[x^2-y^2,y^2-z^2,x^2-z^2])
  @test length(minimal_generating_set(J2)) == 1
  @test ideal(Q2, minimal_generating_set(J2)) == J2
  J3 = ideal(Q3,[x^2-z^2,y*(y-1)])
  @test length(small_generating_set(J3)) == 1
  @test ideal(Q3,small_generating_set(J3)) == J3

  R, (x, y) = QQ[:x, :y]
  I = ideal(R, x^6-y)
  U = complement_of_point_ideal(R, [1, 1])
  L = Oscar.MPolyQuoLocRing(R, I, U)
  J = ideal(L, [x-1, y-1])^2
  minJ = minimal_generating_set(J)
  @test length(minJ) == 1
  @test ideal(L,minJ) == J
end

@testset "dimensions of localizations at prime ideals" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  A, pr = quo(R, ideal(R, [x^2 + 1]))
  U = complement_of_prime_ideal(modulus(A) + ideal(R, y))
  L, loc = localization(A, U)
  @test dim(L) == 1

  U = complement_of_prime_ideal(modulus(A))
  L, loc = localization(A, U)
  @test dim(L) == 0

  U = complement_of_prime_ideal(modulus(A) + ideal(R, [y, z]))
  L, loc = localization(A, U)
  @test dim(L) == 2
end

@testset "composition of morphisms with coefficient maps" begin
  R, (x, y) = ZZ[:x, :y]
  U = powers_of_element(x)
  L, _ = localization(R, U)

  Pt, t = QQ[:t]
  KK, I = number_field(t^2 + 1)

  O_KK = maximal_order(KK)

  S, (u, v) = KK[:u, :v]
  Sl, _ = localization(S, powers_of_element(u))

  R_to_L = hom(R, L, L.([x, y]))
  R_to_Sl = hom(R, Sl, QQ, [u, v])
  L_to_Sl = hom(L, Sl, R_to_Sl)
  @test Oscar._has_coefficient_map(L_to_Sl)
  compose(R_to_L, L_to_Sl)

  id_L = identity_map(L)
  phi = compose(id_L, L_to_Sl)
  @test phi == L_to_Sl
  psi = compose(L_to_Sl, identity_map(Sl))
  @test psi == L_to_Sl

  prep = hom(R, L, coefficient_ring(R), gens(L))
  id_L = hom(L, L, prep)
  phi = compose(id_L, L_to_Sl)
  @test phi == L_to_Sl
  psi = compose(L_to_Sl, identity_map(Sl))
  @test psi == L_to_Sl

  prep = hom(R, L, identity_map(coefficient_ring(R)), gens(L))
  id_L = hom(L, L, prep)
  phi = compose(id_L, L_to_Sl)
  @test phi == L_to_Sl
  psi = compose(L_to_Sl, identity_map(Sl))
  @test psi == L_to_Sl

  R_to_Sl2 = hom(R, Sl, Sl, [u, v])
  L_to_Sl2 = hom(L, Sl, R_to_Sl2)
  @test Oscar._has_coefficient_map(L_to_Sl2)

  phi = compose(identity_map(L), L_to_Sl2)
  @test phi == L_to_Sl2
  psi = compose(L_to_Sl2, identity_map(Sl))
  @test psi == L_to_Sl2

  R_to_Sl3 = hom(R, Sl, O_KK, [u, v])
  L_to_Sl3 = hom(L, Sl, R_to_Sl3)
  @test Oscar._has_coefficient_map(L_to_Sl3)

  phi = compose(identity_map(L), L_to_Sl3)
  @test phi == L_to_Sl3
  psi = compose(L_to_Sl3, identity_map(Sl))
  @test psi == L_to_Sl3


  FF = GF(7)
  T, (a, b) = FF[:a, :b]
  Tl, _ = localization(T, powers_of_element(a))

  R_to_Tl = hom(R, Tl, FF, [a, b])
  L_to_Tl = hom(L, Tl, R_to_Tl)
  @test Oscar._has_coefficient_map(L_to_Tl)

  dirty_reduction = MapFromFunc(KK, GF(7), x->begin c = coordinates(x); FF(numerator(c[1]))/FF(denominator(c[1])) end)

  S_to_Tl = hom(S, Tl, dirty_reduction, [a, b])
  Sl_to_Tl = hom(Sl, Tl, S_to_Tl)

  id_Sl = identity_map(Sl)
  phi = compose(identity_map(Sl), Sl_to_Tl)
  @test phi == Sl_to_Tl
  psi = compose(Sl_to_Tl, identity_map(Tl))
  @test psi == Sl_to_Tl


  L_to_Tl = compose(L_to_Sl, Sl_to_Tl)
  @test L_to_Tl(one(L)) == one(Tl)

  phi = compose(identity_map(L), L_to_Tl)
  @test phi == L_to_Tl
  psi = compose(L_to_Tl, identity_map(Tl))
  @test psi == L_to_Tl

  L_to_Tl2 = compose(L_to_Sl2, Sl_to_Tl)
  @test L_to_Tl(one(L)) == one(Tl)

  phi = compose(identity_map(L), L_to_Tl2)
  @test phi == L_to_Tl2
  psi = compose(L_to_Tl2, identity_map(Tl))
  @test psi == L_to_Tl2

  L_to_Tl3 = compose(L_to_Sl3, Sl_to_Tl)
  @test L_to_Tl(one(L)) == one(Tl)

  phi = compose(identity_map(L), L_to_Tl3)
  @test phi == L_to_Tl3
  psi = compose(L_to_Tl3, identity_map(Tl))
  @test psi == L_to_Tl3


  dirty_reduction2 = MapFromFunc(KK, T, x->begin c = coordinates(x); T(numerator(c[1]))/T(denominator(c[1])) end)

  S_to_Tl2 = hom(S, Tl, dirty_reduction2, [a, b])
  Sl_to_Tl2 = hom(Sl, Tl, S_to_Tl2)

  L_to_Tl = compose(L_to_Sl, Sl_to_Tl2)
  @test L_to_Tl(one(L)) == one(Tl)

  L_to_Tl2 = compose(L_to_Sl2, Sl_to_Tl2)
  @test L_to_Tl(one(L)) == one(Tl)

  L_to_Tl3 = compose(L_to_Sl3, Sl_to_Tl2)
  @test L_to_Tl(one(L)) == one(Tl)


  dirty_reduction3 = MapFromFunc(KK, Tl, x->begin c = coordinates(x); Tl(numerator(c[1]))/Tl(denominator(c[1])) end)

  S_to_Tl3 = hom(S, Tl, dirty_reduction3, [a, b])
  Sl_to_Tl3 = hom(Sl, Tl, S_to_Tl3)

  L_to_Tl = compose(L_to_Sl, Sl_to_Tl3)
  @test L_to_Tl(one(L)) == one(Tl)

  L_to_Tl2 = compose(L_to_Sl2, Sl_to_Tl3)
  @test L_to_Tl(one(L)) == one(Tl)

  L_to_Tl3 = compose(L_to_Sl3, Sl_to_Tl3)
  @test L_to_Tl(one(L)) == one(Tl)
end

@testset "modulus - MPAnyQuoRing MPAnyNonQuoRing" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, [x+y+z])
  A, _ = quo(R,I)
  @test modulus(R) == ideal(R,[zero(R)])
  @test modulus(A) == I
  U= Oscar.MPolyComplementOfKPointIdeal(R,[0,0,0])
  Rl,_ = localization(R,U)
  Il = Rl(I)
  Al, _ = quo(Rl, Il)
  @test modulus(Rl) == ideal(Rl,[zero(Rl)])
  @test modulus(Al) == Il
  U2= Oscar.MPolyComplementOfPrimeIdeal(ideal(R,[x^2+1,y-x,z]))
  Rl2,_ = localization(R,U2)
  Il2 = Rl2(I)
  Al2,_ = quo(Rl2,Il2)
  @test modulus(Rl2) == ideal(Rl2,[zero(Rl2)])
  @test modulus(Al2) == Il2
  U3= Oscar.MPolyPowersOfElement(x+y)
  Rl3,_ = localization(R,U3)
  Il3 = Rl3(I)
  Al3,_ = quo(Rl3,Il3)
  @test modulus(Rl3) == ideal(Rl3,[zero(Rl3)])
  @test modulus(Al3) == Il3
end

@testset "MPolyAnyIdeal printing" begin
  R, x = polynomial_ring(QQ, 100, :x)
  @test Oscar._get_generators_string_one_line(ideal(R, [x[1]])) == "(x1)"
  @test Oscar._get_generators_string_one_line(ideal(R, [x[1], x[2]])) == "(x1, x2)"
  @test Oscar._get_generators_string_one_line(ideal(R, x)) == "with 100 generators"
  @test Oscar._get_generators_string_one_line(ideal(R, [sum(x)])) == "with 1 generator"
end

@testset "monomial_basis MPolyQuoLocRing" begin
  R, (x,y) = QQ["x", "y"]
  L, _ = localization(R, complement_of_point_ideal(R, [0,0]))
  Q1, _ = quo(L, ideal(L, L(x+1)))
  @test isempty(monomial_basis(Q1))
  Q2, _ = quo(L, ideal(L, L(x^2)))  
  @test_throws AbstractAlgebra.InfiniteDimensionError(check_available = true) monomial_basis(Q2)
  Q3, _ = quo(L, ideal(L, L.([x^2, y^3])))  
  @test monomial_basis(Q3) == L.([x*y^2, y^2, x*y, y, x, 1])
  Q4,_ = quo(L, ideal(L, [x^2-x, y^2-2*y]))   
  @test monomial_basis(Q4) == [L(1)]          #test for difference in localized and non-localized case
  @test Oscar._monomial_basis(L, ideal(L, [x^3*(x-1), y*(y-1)*(y-2)])) == L.([x^2, x, 1])
  L1, _ = localization(R, complement_of_point_ideal(R, [1,2]))
  @test Oscar._monomial_basis(L1, ideal(L1, [(x-1)^2, (y-2)^2])) == L1.([x*y, y, x, 1])  
  @test isempty(Oscar._monomial_basis(L1, ideal(L1, L1.([x, y]))))
end

@testset "dimensions" begin
  # to address issue #2721
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, x)
  
  Q, pr = quo(R, I)
  
  W, loc = localization(Q, complement_of_prime_ideal(I))
  J = ideal(W, x)
  @test dim(J) == 0
  K = ideal(W, x-1)
  @test dim(K) == -inf
  
  WW, _ = quo(W, K)
  @test dim(WW) == -inf
  @test dim(underlying_quotient(WW)) == -inf
end

@testset "simplification of subquotients" begin
  R,(x,y) = polynomial_ring(GF(3),2)
  I = ideal(x)
  A, _ = quo(R, I)
  U = complement_of_point_ideal(R, [0, 0])
  L, _ = localization(A, U)
  AA, iso, iso_inv = simplify(L)
  @test all(x==iso(iso_inv(x)) for x in gens(AA))
  @test all(x==iso_inv(iso(x)) for x in gens(L))
end

