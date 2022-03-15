@testset "mpoly-localizations" begin
  R, var = ZZ["x", "y"]
  x = var[1]
  y = var[2] 
  f = x^2 + y^2 -1
  m = ideal(R, [x, y])
  I = ideal(R, f)
  S = MPolyComplementOfPrimeIdeal(I)
  V = Localization(S)
  T = MPolyComplementOfKPointIdeal(R, [ZZ(1), ZZ(0)])
  W = Localization(T)
  
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
  V = Localization(S)
  @test base_ring(V) == R
  @test inverted_set(V) == S

  T = MPolyComplementOfKPointIdeal(R, [k(p), k(q)])
  @test ambient_ring(T) == R
  @test typeof(point_coordinates(T)) == Vector{elem_type(k)}
  @test x in T
  @test !(x-p in T)
  W = Localization(T)

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
  lbpa = groebner_basis(J)
  @test MonomialOrdering(R, ordering(lbpa)) == negdegrevlex(gens(base_ring(W)))
  @test oscar_gens(lbpa)[1] == W(1)

  J = ideal(W, [f, y-q])
  lbpa = groebner_basis(J)
  @test MonomialOrdering(R, ordering(lbpa)) == negdegrevlex(gens(base_ring(W)))
  @test oscar_gens(lbpa) == W.([x-p, y-q])

  @test I_loc + J isa MPolyLocalizedIdeal
  @test I_loc * J isa MPolyLocalizedIdeal

  @test reduce(W(x)//W(y-q+1), lbpa) == W(p)//W(y-q+1)

  K = ideal(W, f)
  @test f*(x-p+4) in K
  @test !(f+2 in K)
  @test f*(x-p+9) in I_loc

  K = ideal(V, [f*x, f*y])
  lbpa = groebner_basis(K)
  reduce(V(x), lbpa)
  K = ideal(V, [x, y])
  lbpa = groebner_basis(K)
  reduce(V(x), lbpa)
  
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
  W = Localization(S)
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

  U = Localization(R, f)
  V = Localization(U, x)
  W = Localization(V, 5*f)

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
# test_Ring_interface_recursive(Localization(S))
# test_Ring_interface_recursive(Localization(T))
# test_Ring_interface_recursive(Localization(U))
  
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

  test_Ring_interface_recursive(Localization(S))
  test_Ring_interface_recursive(Localization(T))
  test_Ring_interface_recursive(Localization(U))

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
# test_Ring_interface_recursive(Localization(S))
# test_Ring_interface_recursive(Localization(T))
# test_Ring_interface_recursive(Localization(U))
end
