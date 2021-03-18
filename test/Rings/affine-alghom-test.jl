@testset "algebra homomorphisms" begin
  r, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  s, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"])
  S = quo(s, ideal(s, [c-b^3]))[1]
  V = S.([2*a+b^6, 7*b-a^2, c^2])
  f = hom(r, S, V)
  K = kernel(f)
  R = quo(r, K)[1]
  phi = hom(R, S, V)
  psi = inverse(phi)

  @test issurjective(f) == true
  @test isinjective(f) == false
  @test isbijective(f) == false
  @test isfinite(f) == true
  @test issurjective(phi) == true
  @test isinjective(phi) == true
  @test isbijective(phi) == true
  @test isfinite(phi) == true
end

