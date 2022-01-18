@testset "mpolyquo-localizations" begin
  #kk = GF(101)
  kk = QQ
  R, (x,y) = kk["x", "y"]

  f = x^2 + y^2-2
  S = MPolyPowersOfElement(x-1)
  T = MPolyComplementOfKPointIdeal(R, [1,1])
  V = MPolyComplementOfPrimeIdeal(ideal(R, f))
  @test S ⊂ V
  @test !(V ⊂ S)
  @test !(T ⊂ V)
  @test (V ⊂ T)
  @test !(MPolyComplementOfPrimeIdeal(ideal(R, f-1)) ⊂ T)
  @test S ⊂ MPolyComplementOfPrimeIdeal(ideal(R, f-1))
  @test !(MPolyPowersOfElement(f) ⊂ V)
  @test MPolyPowersOfElement(x-1) ⊂ MPolyComplementOfKPointIdeal(R, [0,0])
  @test MPolyPowersOfElement(x-1) * MPolyComplementOfKPointIdeal(R, [0,0]) ⊂ MPolyComplementOfKPointIdeal(R, [0,0])
  @test !(MPolyPowersOfElement(x) ⊂ MPolyComplementOfKPointIdeal(R, [0,0]))
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
  W = Localization(U)
  Localization(W, S)
  @test base_ring(W) == R
  @test inverted_set(W) == U
  L = quo(W, ideal(W, f))
  @test issubset(modulus(L), saturated_ideal(localized_modulus(L)))
  @test gens(L) == L.(gens(R))
  @test x//(y*(x-1)) in L

  I = ideal(L, (x-1)*(y-1))
  @test one(W) in I
  @test isunit(L(y-1))

  h = x^4+23*x*y^3-15
  Q, _ = quo(R, f)
  T = MPolyPowersOfElement(h^3)
  W = Localization(Q, T)
  @test x//(h+3*f) in W
  @test W(x//(h+3*f)) == W(x//h)
  g = W.([rand(R, 0:5, 0:2, 0:1) for i in 1:10])
  @test isone(W(h)*inv(W(h)))
  g = [W(8*x*h+4), W(x, 45*h^2)]
  @test one(localized_ring(W)) in ideal(W, g)

  @test one(W) == dot(write_as_linear_combination(one(W), g), g)


  h = (x+5)*(x^2+10*y)+(y-7)*(y^2-3*x)
  Q, _ = quo(R, h)
  T = MPolyComplementOfKPointIdeal(R, [-5, 7])
  W = Localization(Q, T)
  @test x//(y) in W
  @test x//(y+h) in W
  g = [W(h + (x+5) - 9, y+24*x^3-8)]
  (d, a) = bring_to_common_denominator([one(W), g[1]])
  @test W(a[1], d) == one(W)
  @test W(a[2]*lifted_numerator(g[1]), d) == g[1]
  @test (one(localized_ring(W)) in ideal(W, g))
  @test one(W) == dot(write_as_linear_combination(one(W), g), g) 
end
