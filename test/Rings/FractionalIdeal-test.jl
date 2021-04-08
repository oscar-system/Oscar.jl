@testset "FractionalIdeal" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

  a = fractional_ideal(ideal(R, [x, y]), z)

  @test isa(a, Oscar.FractionalIdeal)

  @test a + a == a

  @test a*a == fractional_ideal(ideal(R, [x^2, x*y, y^2]), z^2)

  @test a == fractional_ideal(numerator(a), denominator(a))

  b = fractional_ideal(ideal(R, [y, z]), x)

  @test a + b == fractional_ideal(ideal(R, [x^2, x*y, z*y, z^2]), x*z)
end
