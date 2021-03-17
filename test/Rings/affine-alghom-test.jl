@testset "algebra homomorphisms" begin
  r, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"])
  s, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"])
  S = quo(s, ideal(s, [c-b^3]))[1]
  V = S.([2*a+b^6, 7*b-a^2, c^2, zero(s)])
  phi = hom(r, S, V)
  K = kernel(phi)
  R = quo(r, K)[1]
  psi = hom(R, S, V)
  
  @test issurjective(phi) == true
  @test isinjective(phi) == false
  @test isbijective(phi) == false
  @test issurjective(psi) == true
  @test isinjective(psi) == true
  @test isbijective(psi) == true
end
