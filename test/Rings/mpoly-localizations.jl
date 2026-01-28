const rng = Oscar.get_seeded_rng()

@testset "mpoly-localizations" begin
  R, variab = ZZ[:x, :y]
  x = variab[1]
  y = variab[2] 
  f = x^2 + y^2 -1
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = Oscar.MPolyComplementOfPrimeIdeal(I)
  V, _ = localization(S)
  
  k = QQ
  R, variab = k[:x, :y]
  x = variab[1]
  y = variab[2] 
  p = 123
  q = -49
  f = x^2 + y^2 - p^2 - q^2
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = Oscar.MPolyComplementOfPrimeIdeal(I)
  @test ring(S) == R
  @test !( f in S )
  @test x in S
  V, _ = localization(S)
  @test base_ring(V) == R
  @test inverted_set(V) == S

  T = Oscar.MPolyComplementOfKPointIdeal(R, [k(p), k(q)])
  @test ring(T) == R
  @test typeof(point_coordinates(T)) == Vector{elem_type(k)}
  @test x in T
  @test !(x-p in T)
  W, _ = localization(T)

  a = W(R(1))
  b = W(2)
  c = W(y//x)
  @test a*(b+c) == a*b + a*c
  @test a*(b+c)^2 != a*(b^2) + a*c^2
  b = W((x-19)//(y-2))
  W(1)/b
  c = W((f+3)//(y-9))
  @test c/b == 1/(b/c)

  I_loc = ideal(W, gens(I))
  @test gens(I_loc) == [W(f)]
  I_loc = ideal(W, gen(I, 1))
  @test gens(I_loc) == [W(f)]
  #@test W(f) in I_loc
  
  J = ideal(W, [W(x-3), W(y-4)])

  J = ideal(W, [f, y-q])
  @test I_loc + J isa Oscar.MPolyLocalizedIdeal
  @test I_loc * J isa Oscar.MPolyLocalizedIdeal

 # @test reduce(W(x)//W(y-q+1), lbpa) == W(p)//W(y-q+1)

  K = ideal(W, f)
  @test f*(x-p+4) in K
  @test !(f+2 in K)
  @test f*(x-p+9) in I_loc

  K = ideal(V, [f*x, f*y])
  K = ideal(V, [x, y])
  
  R, v = ZZ[:x, :y]
  x = v[1]
  y = v[2] 
  f = (x^2 + y^2)^2
  T = Oscar.MPolyPowersOfElement(R, [f])
  S = Oscar.MPolyComplementOfPrimeIdeal(ideal(R, [x, y]))
  U = Oscar.MPolyProductOfMultSets(R, [S, T])
  @test f in U
  @test (f*(x-1) in U)
  @test !(f*x in U)

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  o = degrevlex([x, y])*negdegrevlex([z])
  S, _ = localization(R, o)
  @test z + 1 in inverted_set(S)
  @test !(x + 1 in inverted_set(S))
  I = ideal(S, [x + y + z, x^2 + y^2 + z^3])
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)

  R, (x, y) = QQ[:x, :y]
  U = Oscar.MPolyPowersOfElement(R, [x])
  L = Oscar.MPolyLocRing(R, U)
  I = ideal(L, [y*(y-x)])
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)
  J = ideal(L, [x*(y-x)])
  @test y in quotient(I, J)
  @test I:J == quotient(I, J)

  R, (x, y) = QQ[:x, :y]
  f = x^2 + y^3- 2
  I = ideal(R, [f])
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)
  U = Oscar.MPolyComplementOfPrimeIdeal(I)
  L = Oscar.MPolyLocRing(R, U)
  I = ideal(L, [f^4])
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)
  J = ideal(L, [f^2])
  @test f^2 in quotient(I, J)
  @test !(f in quotient(I, J))
  @test I:J == quotient(I, J)
  
  R, (x, y) = QQ[:x, :y]
  f = x^2 + y^2- 2
  I = ideal(R, [f])
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)
  U = Oscar.MPolyComplementOfKPointIdeal(R, [1, 1])
  L = Oscar.MPolyLocRing(R, U)
  I = ideal(L, [y^2*f^4])
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)
  J = ideal(L, [x*f^2])
  @test f^2 in quotient(I, J)
  @test !(f in quotient(I, J))
  @test I:J == quotient(I, J)
end

@testset "mpoly-localizations PowersOfElements" begin
  R, variab = ZZ[:x, :y]
  x = variab[1]
  y = variab[2] 
  f = x^2 + y^2 -1
  S = Oscar.MPolyPowersOfElement(R, [x, y, f])
  @test f in S
  # 5 is not a unit in R
  @test !(5*f in S)

  R, variabs = QQ[:x, :y]
  x = variabs[1]
  y = variabs[2]
  f = x^2 + y^4-120
  S = Oscar.MPolyPowersOfElement(R, [x, y, f])
  @test f in S
  # Now 5 is a unit in R
  @test (5*f in S)
  @test x^3*y in S
  @test !(x^19*(x+y) in S)
  W, _ = localization(S)
  @test x*y^2 in S
  @test !(x*(x+y) in S)
  @test W(1//x) == 1/W(x)
  I = ideal(W, x*y*(x+y))
  small_generating_set(I; algorithm = :simple)
  small_generating_set(I; algorithm = :with_saturation)
  @test x+y in I
  @test x+y in saturated_ideal(I)
  @test !(x+2*y in I)
  J = I + ideal(W, f)
  @test W(1) in J
end

@testset "mpoly-localization homomorphisms" begin
  R, variab = ZZ[:x, :y]
  x = variab[1]
  y = variab[2] 
  f = x^2 + y^2 -1
  S = Oscar.MPolyPowersOfElement(R, [x, y, f])
  @test f in S
  @test !(5*f in S)

  U, _ = localization(R, f)
  V, _ = localization(U, x)
  W, _ = localization(V, 5*f)

  phi = Oscar.MPolyLocalizedRingHom(R, W, [W(1//x), W(y//f)])
  @test phi(x*f) == phi(domain(phi)(x*f))
  @test phi(ZZ(9)) == W(ZZ(9))

  psi = Oscar.MPolyLocalizedRingHom(U, R, [R(0), R(0)])
  @test psi(U(-1//f))==one(codomain(psi))
end

@testset "Ring interface for localized polynomial rings" begin
# kk = QQ
# R, v = kk[:x, :y]
# x = v[1]
# y = v[2] 
# f = x^2+y^2-1
# I = ideal(R, f)
#
# d = Vector{elem_type(R)}()
# d = [rand(R, 1:3, 0:4, 1:10)::elem_type(R) for i in 0:(abs(rand(Int))%3+1)]
# S = Oscar.MPolyPowersOfElement(R, d)
# T = Oscar.MPolyComplementOfKPointIdeal(R, [kk(125), kk(-45)])
# U = Oscar.MPolyComplementOfPrimeIdeal(I)
#
# ConformanceTests.test_Ring_interface_recursive(localization(S)[1])
# ConformanceTests.test_Ring_interface_recursive(localization(T)[1])
# ConformanceTests.test_Ring_interface_recursive(localization(U)[1])

# should be unnecessary: https://github.com/oscar-system/Oscar.jl/pull/1459#issuecomment-1230185617
#  AbstractAlgebra.promote_rule(::Type{fpMPolyRingElem}, ::Type{ZZRingElem}) = fpMPolyRingElem
#  AbstractAlgebra.promote_rule(::Type{fpFieldElem}, ::Type{ZZRingElem}) = fpFieldElem

  kk = GF(7) 
  R, v = kk[:x, :y]
  x = v[1]
  y = v[2] 
  f = x^2+y^2-1
  I = ideal(R, f)

  d = Vector{elem_type(R)}()
  for i in 0:(abs(rand(rng, Int))%3+1)
    f = rand(rng, R, 1:3, 0:3, 1:10)::elem_type(R)
    iszero(f) || push!(d, f)
  end
  S = Oscar.MPolyPowersOfElement(R, d)
  T = Oscar.MPolyComplementOfKPointIdeal(R, [kk(125), kk(-45)])
  U = Oscar.MPolyComplementOfPrimeIdeal(I)

  ConformanceTests.test_Ring_interface_recursive(localization(S)[1])
  ConformanceTests.test_Ring_interface_recursive(localization(T)[1])
  ConformanceTests.test_Ring_interface_recursive(localization(U)[1])

# kk = ZZ
# R, v = kk[:x, :y]
# x = v[1]
# y = v[2] 
# f = x^2+y^2-1
# I = ideal(R, f)
#
# d = Vector{elem_type(R)}()
# d = [rand(R, 1:3, 0:4, 1:10)::elem_type(R) for i in 0:(abs(rand(Int))%3+1)]
# S = Oscar.MPolyPowersOfElement(R, d)
# U = Oscar.MPolyComplementOfPrimeIdeal(I)
#
# ConformanceTests.test_Ring_interface_recursive(localization(S)[1])
# ConformanceTests.test_Ring_interface_recursive(localization(T)[1])
# ConformanceTests.test_Ring_interface_recursive(localization(U)[1])
end

@testset "localization_at_orderings_1" begin
  R, (x,y) = QQ[:x, :y]
  o = degrevlex([x])*negdegrevlex([y])
  U = Oscar.MPolyLeadingMonOne(R, o)
  @test y-1 in U
  @test !(x in U)
  L, _ = localization(R, U)
  I = ideal(L, [x^2, y*(y-1)])
  @test !(x in I)
  @test x^2 in I
  @test y in I
  @test dot(coordinates(y, I), gens(I)) == y
end

@testset "localization_at_orderings_2" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  o = degrevlex([x, y])*negdegrevlex([z])
  S, _ = localization(R, o)
  @test z + 1 in inverted_set(S)
  @test !(x + 1 in inverted_set(S))
  I = ideal(S, [x + y + z, (z+1)*(x^2 + y^2 + z^3)])
  f = (x^2 + y^2 + z^3) + 5*(x + y + z)
  @test f in I
  @test coordinates(f, I) == S[5 1//(z+1)]
end

@testset "localizations at k-points" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  p = [-5, 8, 1//2]
  U = Oscar.MPolyComplementOfKPointIdeal(R, p)
  I = ideal(R, [x*(x+5), (y-8)*y-z*(x+5)])
  L, _ = localization(R, U)
  LI = L(I)
  @test x+5 in LI
  @test y-8 in LI
  @test dot(coordinates(y-8, LI), gens(LI)) == L(y-8)
  @test dot(coordinates(x+5, LI), gens(LI)) == L(x+5)

  W, _ = quo(L, LI)
  J = ideal(W, [(z-1//2)^4*y])
  @test (z-1//2)^5 in J
  @test !((z-1//2)^5 in LI)
  @test coordinates((z-1//2)^5, J) == matrix(W, 1, 1, [W(z-1//2, y)])
end

@testset "successive localizations" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  p = [0,0,0]
  U = Oscar.MPolyComplementOfKPointIdeal(R, p)
  I = ideal(R, [x*(y-1)-z*(x-2), y*x])
  L, _ = localization(R, U)
  LI = L(I)
  W, _ = quo(L, LI)
  S = Oscar.MPolyPowersOfElement(R, [y])
  RS, _ = localization(R, S)
  RSI = RS(I)
  saturated_ideal(RSI, with_generator_transition=true)
  J = L(Oscar.pre_saturated_ideal(RSI))
  z in J
  W, _ = quo(L, LI)
  S = Oscar.MPolyPowersOfElement(R, [y])
  WS, _ = localization(W, S)
  @test !iszero(W(z))
  @test iszero(WS(z))
  LS, _ = localization(L, S)
  LSI = LS(LI)
  @test dot(coordinates(z, LSI), gens(LSI)) == LS(z)
  @test !(z in Oscar.pre_saturated_ideal(LSI))
  @test z in LSI
  @test !(z in Oscar.pre_saturated_ideal(LSI)) # caching is not supposed to happen, because of special routing.
end

@testset "saturated_ideals" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, [x, y^2+1])
  U = Oscar.MPolyComplementOfPrimeIdeal(I)
  L = Oscar.MPolyLocRing(R, U)
  J = ideal(L,[y*(x^2+(y^2+1)^2)])
  J_sat = ideal(R,[(x^2+(y^2+1)^2)])
  @test saturated_ideal(J) == J_sat
  @test_throws ErrorException("computation of the transition matrix for the generators is not supposed to happen for localizations at complements of prime ideals") saturated_ideal(L(J_sat); with_generator_transition=true)
  JJ = ideal(R,[y*(x^2+(y^2+1)^2)])
  @test saturated_ideal(JJ) == JJ
end

@testset "zero divisors" begin
  R, (x, y) = QQ[:x, :y]
  @test !is_zero_divisor(x)
  @test !is_zero_divisor(R(5))
  @test is_zero_divisor(zero(x))
  I = ideal(R, x*y)
  A, _ = quo(R, I)
  @test !is_zero_divisor(A(x+y))
  @test !is_zero_divisor(A(5))
  @test is_zero_divisor(A(x))
  L, _ = localization(R, x+y)
  @test !is_zero_divisor(L(5))
  @test !is_zero_divisor(L(x))
  @test is_zero_divisor(zero(L))
  J = ideal(R, (x+y)*x*y)
  A2, _ = quo(R, J)
  W, _ = localization(A2, A2(x+y))
  @test !is_zero_divisor(W(x-y))
  @test is_zero_divisor(W((x+y)^7*x))

  Z4, _ = quo(ZZ, 4)
  Z4x, (x, y) = Z4[:x, :y]
  f = 2*x
  @test is_zero_divisor(Z4(2))
  @test is_zero_divisor(Z4x(2))
  @test is_zero_divisor(f)
  # At the moment, the Singular backend complains when the modulus of the coefficient 
  # ring is not a prime. So we leave the following tests out for the time being. 
#  L, _ = localization(Z4x, Z4x(x))
#  @test !is_zero_divisor(L(5))
#  @test is_zero_divisor(L(2))
#  @test !is_zero_divisor(L(x))
#  @test is_zero_divisor(2*L(x))
#  @test is_zero_divisor(zero(L))
#  I = ideal(Z4x, x*y)
#  A, _ = quo(Z4x, I)
#  @test !is_zero_divisor(A(x+y))
#  @test is_zero_divisor(2*A(x+y))
#  @test !is_zero_divisor(A(5))
#  @test is_zero_divisor(A(2))
#  @test is_zero_divisor(A(x))
#  W, _ = localization(A, A(x-y))
#  @test !is_zero_divisor(W(x+y))
#  @test is_zero_divisor(2*W(x+y))
#  @test !is_zero_divisor(W(5))
#  @test is_zero_divisor((x-y)*W(2))
#  @test is_zero_divisor((x-y)*W(x))
end

@testset "saturated ideals II" begin
  R, (x,y) = QQ[:x, :y]
  L, _ = localization(R, powers_of_element(x))
  I = ideal(R, [x^7*(y-1), (x^5)*(x-1)])
  J = L(I)
  @test (x-1) in J
  c1 = coordinates(x-1, J)
  @test x-1 in saturated_ideal(J)
  c2 = coordinates(y-1, J)
  @test y-1 == c2[1]*gen(J, 1) + c2[2]*gen(J, 2)
  c3 = coordinates(x-y, J)
  @test x-y == c3[1]*gen(J, 1) + c3[2]*gen(J, 2)
  JJ = L(I)
  @test x-1 in saturated_ideal(JJ, strategy=:single_saturation)
  c2 = coordinates(y-1, JJ)
  @test y-1 == c2[1]*gen(JJ, 1) + c2[2]*gen(JJ, 2)
end

@testset "printing" begin
  R, (x,y) = ZZ[:x, :y]
  U = powers_of_element(x)
  L, m = localization(R, U)
  @test sprint(show, L(y, x)) == "y/x"
  @test sprint(show, L(y)) == "y" # Denominator is one; omit printing it.
end

@testset "fractions and divexact: Issue #1885" begin
  R, (x, y) = ZZ[:x, :y]
  U = powers_of_element(y)
  L, = localization(R, U)

  @test L(x)/L(y) == L(x//y)
  @test_throws ErrorException L(y)/L(x-1)
end

@testset "minimal generating sets" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, [x-1, y-2])^2
  J = ideal(R, [x, y])^2
  L1, phi1 = localization(R, complement_of_point_ideal(R, [1, 2]))
  L2, phi2 = localization(R, powers_of_element(x))

  @test length(minimal_generating_set(phi1(I))) == 3
  @test length(small_generating_set(phi1(I))) == 3
  @test length(small_generating_set(phi2(I))) == 3

  @test length(minimal_generating_set(phi1(J))) == 1
  @test isone(first(minimal_generating_set(phi1(J))))

  small_gens = small_generating_set(phi2(J))
  @test length(small_gens) == 3 # saturated_ideal is not called here.
  #@test isone(first(small_gens))
end

@testset "dimensions of localizations at prime ideals" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  f = x^2 + x*y^2 - z^3
  P = ideal(R, f)
  U = complement_of_prime_ideal(P)
  L, loc = localization(R, U)
  @test dim(L) == 1

  U = complement_of_prime_ideal(ideal(R, [x^2 + 1, y]))
  L, loc = localization(R, U)
  @test dim(L) == 2

  U = complement_of_prime_ideal(ideal(R, [x^2 + 1, y, z]))
  L, loc = localization(R, U)
  @test dim(L) == 3
end

@testset "dimensions of localized ideals" begin
  R, (x, y) = QQ[:x, :y]
  U = complement_of_point_ideal(R, [0, 1])

  L, loc_map = localization(R, U)

  I = ideal(R, [y*x, y*(y-1)])
  @test dim(I) == 1
  J = loc_map(I)
  @test dim(J) == 0
end

@testset "radical computation" begin
  R, (x, y) = QQ[:x, :y]
  #U = complement_of_prime_ideal(ideal(R, y))
  U = powers_of_element(x)
  J = ideal(R, x-y^2+ 1)
  W, _ = quo(R, J)

  L, loc_map = localization(W, U)

  I = ideal(L, [y^2, 0]) # Dealing with zero generators lead to a bug once, so this is a real test.
  @test radical(I) == ideal(L, [y])
end

@testset "mpoly-loc constructors" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  m = ideal(R, [y - 1, x - 2, z - 3])
  Q = localization(R, m)
  I = ideal(Q, [x - 2, (y - 1)^2*z])
  a = Q(x//(1 + x))

  @test gens(I) == Q.([x - 2, y^2*z - 2*y*z + z])
  @test fraction(a) == x // (1 + x)
end

@testset "mpoly-loc operations" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  m = ideal(R, [y - 1, x - 2, z - 3])
  Q = localization(R, m)
  I = ideal(Q, [x - 2, (y - 1)^2*z])
  J = ideal(Q, [x - 2, y - 1])
  a = Q(x//(1 + x))
  b = Q((x + y)//(3 + x - y))

  K1 = I+J
  K2 = I*J
  K3 = J^3
  c1 = a+b
  c2 = a*b

  @test iszero(a-a)
  @test a*deepcopy(a) == a^2
  @test det(matrix(Q, [1 x; y z])) == Q(z-x*y)

  @test gens(K1) == Q.([x - 2, y^2*z - 2*y*z + z, x - 2, y - 1])
  @test gens(K2) == Q.([x^2 - 4*x + 4, x*y - x - 2*y + 2, x*y^2*z - 2*x*y*z + x*z - 2*y^2*z + 4*y*z - 2*z, y^3*z - 3*y^2*z + 3*y*z - z])
  @test gens(K3) == Q.([x^3 - 6*x^2 + 12*x - 8, x^2*y - x^2 - 4*x*y + 4*x + 4*y - 4, x^2*y - x^2 - 4*x*y + 4*x + 4*y - 4, x*y^2 - 2*x*y + x - 2*y^2 + 4*y - 2, x^2*y - x^2 - 4*x*y + 4*x + 4*y - 4, x*y^2 - 2*x*y + x - 2*y^2 + 4*y - 2, x*y^2 - 2*x*y + x - 2*y^2 + 4*y - 2, y^3 - 3*y^2 + 3*y - 1])
  @test fraction(c1) == (x*(3 + x - y) + (1 + x)*(x + y))//((1 + x)*(3 + x - y))
  @test fraction(c2) == (x*(x + y))//((1 + x)*(3 + x - y))

  (x, y, z) = map(Q, (x, y, z))
  @test inv(x) + inv(y) == (x + y)*inv(x*y)
  @test_throws Exception x*inv(x - 2)
end

@testset "dimensions" begin
  # to address issue #2721
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, x)
  L, loc = localization(R, complement_of_prime_ideal(I))
  J = ideal(L, x)
  @test dim(J) == 0
  @test dim(L) == 1
  K = ideal(L, y)
  @test dim(K) == -inf
end

@testset "issue 5316" begin
  R, (x, y) = QQ[:x, :y]
  f = x+y
  U = powers_of_element(f)
  V = Oscar.MPolyLeadingMonOne(neglex(gens(R)))
  T = U*V
  L, _ = localization(R, T)
  I = ideal(L, x)
  @test coordinates(x, I)[1] == one(L)
end

