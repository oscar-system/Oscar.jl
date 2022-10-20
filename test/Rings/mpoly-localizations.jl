@testset "mpoly-localizations" begin
  R, var = ZZ["x", "y"]
  x = var[1]
  y = var[2] 
  f = x^2 + y^2 -1
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = MPolyComplementOfPrimeIdeal(I)
  V, _ = Localization(S)
  T = MPolyComplementOfKPointIdeal(R, [ZZ(1), ZZ(0)])
  W, _ = Localization(T)
  
  k = QQ
  R, var = k["x", "y"]
  x = var[1]
  y = var[2] 
  p = 123
  q = -49
  f = x^2 + y^2 - p^2 - q^2
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = MPolyComplementOfPrimeIdeal(I)
  @test ambient_ring(S) == R
  @test !( f in S )
  @test x in S
  V, _ = Localization(S)
  @test base_ring(V) == R
  @test inverted_set(V) == S

  T = MPolyComplementOfKPointIdeal(R, [k(p), k(q)])
  @test ambient_ring(T) == R
  @test typeof(point_coordinates(T)) == Vector{elem_type(k)}
  @test x in T
  @test !(x-p in T)
  W, _ = Localization(T)

  a = W(R(1))
  b = W(2)
  c = W(y//x)
  @test a*(b+c) == a*b + a*c
  @test a*(b+c)^2 != a*(b^2) + a*c^2
  b = W((x-19)//(y-2))
  W(1)//b
  c = W((f+3)//(y-9))
  @test c//b == 1//(b//c)

  I_loc = ideal(W, gens(I))
  @test gens(I_loc) == [W(f)]
  I_loc = ideal(W, gens(I)[1])
  @test gens(I_loc) == [W(f)]
  #@test W(f) in I_loc
  
  J = ideal(W, [W(x-3), W(y-4)])

  J = ideal(W, [f, y-q])
  @test I_loc + J isa MPolyLocalizedIdeal
  @test I_loc * J isa MPolyLocalizedIdeal

 # @test reduce(W(x)//W(y-q+1), lbpa) == W(p)//W(y-q+1)

  K = ideal(W, f)
  @test f*(x-p+4) in K
  @test !(f+2 in K)
  @test f*(x-p+9) in I_loc

  K = ideal(V, [f*x, f*y])
  K = ideal(V, [x, y])
  
  R, v = ZZ["x", "y"]
  x = v[1]
  y = v[2] 
  f = (x^2 + y^2)^2
  S = MPolyComplementOfKPointIdeal(R, [ZZ(0), ZZ(0)])
  T = MPolyPowersOfElement(R, [f])
  U = MPolyProductOfMultSets(R, [S, T])
  @test f in U
  @test (f*(x-1) in U)
  @test !(f*x in U)

  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  o = degrevlex([x, y])*negdegrevlex([z])
  S, _ = Localization(R, o)
  @test z + 1 in inverted_set(S)
  @test !(x + 1 in inverted_set(S))
  I = ideal(S, [x + y + z, x^2 + y^2 + z^3])

  R, (x, y) = QQ["x", "y"]
  U = MPolyPowersOfElement(R, [x])
  L = MPolyLocalizedRing(R, U)
  I = ideal(L, [y*(y-x)])
  J = ideal(L, [x*(y-x)])
  @test y in quotient(I, J)
  @test I:J == quotient(I, J)

  R, (x, y) = QQ["x", "y"]
  f = x^2 + y^3- 2
  I = ideal(R, [f])
  U = MPolyComplementOfPrimeIdeal(I)
  L = MPolyLocalizedRing(R, U)
  I = ideal(L, [f^4])
  J = ideal(L, [f^2])
  @test f^2 in quotient(I, J)
  @test !(f in quotient(I, J))
  @test I:J == quotient(I, J)
  
  R, (x, y) = QQ["x", "y"]
  f = x^2 + y^2- 2
  I = ideal(R, [f])
  U = MPolyComplementOfKPointIdeal(R, [1, 1])
  L = MPolyLocalizedRing(R, U)
  I = ideal(L, [y^2*f^4])
  J = ideal(L, [x*f^2])
  @test f^2 in quotient(I, J)
  @test !(f in quotient(I, J))
  @test I:J == quotient(I, J)
end

@testset "mpoly-localizations PowersOfElements" begin
  R, var = ZZ["x", "y"]
  x = var[1]
  y = var[2] 
  f = x^2 + y^2 -1
  S = MPolyPowersOfElement(R, [x, y, f])
  @test f in S
  # 5 is not a unit in R
  @test !(5*f in S)

  R, vars = QQ["x","y"]
  x = vars[1]
  y = vars[2]
  f = x^2 + y^4-120
  S = MPolyPowersOfElement(R, [x, y, f])
  @test f in S
  # Now 5 is a unit in R
  @test (5*f in S)
  @test x^3*y in S
  @test !(x^19*(x+y) in S)
  W, _ = Localization(S)
  @test x*y^2 in S
  @test !(x*(x+y) in S)
  @test W(1//x) == 1//W(x)
  I = ideal(W, x*y*(x+y))
  @test x+y in I
  @test !(x+2*y in I)
  J = I + ideal(W, f)
  @test W(1) in J
end

@testset "mpoly-localization homomorphisms" begin
  R, var = ZZ["x", "y"]
  x = var[1]
  y = var[2] 
  f = x^2 + y^2 -1
  S = MPolyPowersOfElement(R, [x, y, f])
  @test f in S
  @test !(5*f in S)

  U, _ = Localization(R, f)
  V, _ = Localization(U, x)
  W, _ = Localization(V, 5*f)

  phi = MPolyLocalizedRingHom(R, W, [W(1//x), W(y//f)])
  @test phi(x*f) == phi(domain(phi)(x*f))
  @test phi(ZZ(9)) == W(ZZ(9))

  psi = MPolyLocalizedRingHom(U, R, [R(0), R(0)])
  @test psi(U(-1//f))==one(codomain(psi))
end

function test_elem(W::MPolyLocalizedRing) 
  f = rand(W, 0:3, 0:4, 0:3)
  return f
end

@testset "Ring interface for localized polynomial rings" begin
  include(joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))
#
# kk = QQ
# R, v = kk["x", "y"]
# x = v[1]
# y = v[2] 
# f = x^2+y^2-1
# I = ideal(R, f)
#
# d = Vector{elem_type(R)}()
# d = [rand(R, 1:3, 0:4, 1:10)::elem_type(R) for i in 0:(abs(rand(Int))%3+1)]
# S = MPolyPowersOfElement(R, d)
# T = MPolyComplementOfKPointIdeal(R, [kk(125), kk(-45)])
# U = MPolyComplementOfPrimeIdeal(I)
#
# test_Ring_interface_recursive(Localization(S)[1])
# test_Ring_interface_recursive(Localization(T)[1])
# test_Ring_interface_recursive(Localization(U)[1])
  
  AbstractAlgebra.promote_rule(::Type{gfp_mpoly}, ::Type{fmpz}) = gfp_mpoly
  AbstractAlgebra.promote_rule(::Type{gfp_elem}, ::Type{fmpz}) = gfp_elem

  kk = GF(7) 
  R, v = kk["x", "y"]
  x = v[1]
  y = v[2] 
  f = x^2+y^2-1
  I = ideal(R, f)

  d = Vector{elem_type(R)}()
  for i in 0:(abs(rand(Int))%3+1)
    f = rand(R, 1:3, 0:3, 1:10)::elem_type(R)
    iszero(f) || push!(d, f)
  end
  S = MPolyPowersOfElement(R, d)
  T = MPolyComplementOfKPointIdeal(R, [kk(125), kk(-45)])
  U = MPolyComplementOfPrimeIdeal(I)

  test_Ring_interface_recursive(Localization(S)[1])
  test_Ring_interface_recursive(Localization(T)[1])
  test_Ring_interface_recursive(Localization(U)[1])

# kk = ZZ
# R, v = kk["x", "y"]
# x = v[1]
# y = v[2] 
# f = x^2+y^2-1
# I = ideal(R, f)
#
# d = Vector{elem_type(R)}()
# d = [rand(R, 1:3, 0:4, 1:10)::elem_type(R) for i in 0:(abs(rand(Int))%3+1)]
# S = MPolyPowersOfElement(R, d)
# T = MPolyComplementOfKPointIdeal(R, [kk(125), kk(-45)])
# U = MPolyComplementOfPrimeIdeal(I)
#
# test_Ring_interface_recursive(Localization(S)[1])
# test_Ring_interface_recursive(Localization(T)[1])
# test_Ring_interface_recursive(Localization(U)[1])
end

@testset "localization_at_orderings_1" begin
  R, (x,y) = QQ["x", "y"]
  o = degrevlex([x])*negdegrevlex([y])
  U = MPolyLeadingMonOne(R, o)
  @test y-1 in U
  @test !(x in U)
  L, _ = Localization(R, U)
  I = ideal(L, [x^2, y*(y-1)])
  @test !(x in I)
  @test x^2 in I
  @test y in I
  @test dot(coordinates(y, I), gens(I)) == y
end

@testset "localization_at_orderings_2" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  o = degrevlex([x, y])*negdegrevlex([z])
  S, _ = Localization(R, o)
  @test z + 1 in inverted_set(S)
  @test !(x + 1 in inverted_set(S))
  I = ideal(S, [x + y + z, (z+1)*(x^2 + y^2 + z^3)])
  f = (x^2 + y^2 + z^3) + 5*(x + y + z)
  @test f in I
  @test coordinates(f, I) == S[5 1//(z+1)]
end

@testset "localizations at k-points" begin
  R, (x, y, z) = QQ["x", "y", "z"]
  p = [-5, 8, 1//2]
  U = MPolyComplementOfKPointIdeal(R, p)
  I = ideal(R, [x*(x+5), (y-8)*y-z*(x+5)])
  L, _ = Localization(R, U)
  LI = L(I)
  @test x+5 in LI
  @test y-8 in LI
  @test dot(coordinates(y-8, LI), gens(LI)) == L(y-8)
  @test dot(coordinates(x+5, LI), gens(LI)) == L(x+5)

  W, _ = quo(L, LI)
  J = ideal(W, [(z-1//2)^4*y])
  @test (z-1//2)^5 in J
  @test !((z-1//2)^5 in LI)
  @test coordinates((z-1//2)^5, J) == MatrixSpace(W, 1, 1)([W(z-1//2, y)])
end

@testset "successive localizations" begin
  R, (x, y, z) = QQ["x", "y", "z"]
  p = [0,0,0]
  U = MPolyComplementOfKPointIdeal(R, p)
  I = ideal(R, [x*(y-1)-z*(x-2), y*x])
  L, _ = Localization(R, U)
  LI = L(I)
  W, _ = quo(L, LI)
  S = MPolyPowersOfElement(R, [y])
  RS, _ = Localization(R, S)
  RSI = RS(I)
  saturated_ideal(RSI, with_generator_transition=true)
  J = L(Oscar.pre_saturated_ideal(RSI))
  z in J
  W, _ = quo(L, LI)
  S = MPolyPowersOfElement(R, [y])
  WS, _ = Localization(W, S)
  @test !iszero(W(z))
  @test iszero(WS(z))
  LS, _ = Localization(L, S)
  LSI = LS(LI)
  @test dot(coordinates(z, LSI), gens(LSI)) == LS(z)
  @test !(z in Oscar.pre_saturated_ideal(LSI))
  @test z in LSI
  @test !(z in Oscar.pre_saturated_ideal(LSI)) # caching is not supposed to happen, because of special routing.
end

@testset "zero divisors" begin
  R, (x, y) = QQ["x", "y"]
  @test !is_zero_divisor(x)
  @test is_zero_divisor(zero(x))
  I = ideal(R, x*y)
  A, _ = quo(R, I)
  @test !is_zero_divisor(A(x+y))
  @test is_zero_divisor(A(x))
  L, _ = localization(R, x+y)
  @test !is_zero_divisor(L(x))
  @test is_zero_divisor(zero(L))
  J = ideal(R, (x+y)*x*y)
  A2, _ = quo(R, J)
  W, _ = localization(A2, A2(x+y))
  @test !is_zero_divisor(W(x-y))
  @test is_zero_divisor(W((x+y)^7*x))
end
