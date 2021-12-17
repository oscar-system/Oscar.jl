@testset "mpolyquo-localizations" begin
  #kk = GF(101)
  kk = QQ
  R, (x,y) = kk["x", "y"]

  f = x^2 + y^2-2
  S = MPolyPowersOfElement(x-1)
  T = MPolyComplementOfKPointIdeal(R, [1,1])

  U = S*T
  W = Localization(U)
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
  g = [rand(W, 0:5, 0:2, 0:1) for i in 1:5]
  @test (one(localized_ring(W)) in ideal(W, g)) ? one(W) == dot(write_as_linear_combination(one(W), g), g) : true


  h = (x+5)*(x^2+10*y)+(y-7)*(y^2-3*x)
  Q, _ = quo(R, h)
  T = MPolyComplementOfKPointIdeal(R, [-5, 7])
  W = Localization(Q, T)
  @test x//(y) in W
  @test x//(y+h) in W
  g = [rand(W, 0:5, 0:2, 0:1) for i in 1:5]
  @test (one(localized_ring(W)) in ideal(W, g)) ? one(W) == dot(write_as_linear_combination(one(W), g), g) : true
end
