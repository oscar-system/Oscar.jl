@testset "MPolyQuo" begin
  R, (x,y) = PolynomialRing(QQ, ["x", "y"])

  f = y^2+y+x^2
  C = ideal(R, [f])
  Q, = quo(R, C)

  @test one(Q) == 1
  @test zero(Q) == 0

  I = ideal([x^3 + y^3 - 3, x^5 + y^5 - 5])
  Q, = quo(R, I)
  @test length(Oscar._kbase(Q)) == 12
  b = inv(Q(x))
  @test isone(b*Q(x))

  xx, yy = gens(Q)

  @test xx^2 == xx * xx
  @test xx^0 == one(Q)
end

@testset "MPolyQuo.ideals" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  Q, _ = quo(R, ideal(R, [x*y, x*z]))
  (x, y, z) = map(Q, (x, y, z))

  @test ideal(Q, [x, y, z]) isa Oscar.Ideal

  @test !iszero(ideal(Q, [x, y]))
  @test !iszero(ideal(Q, [y*z]))
  @test iszero(ideal(Q, [x*y, x*z]))
  @test iszero(ideal(Q, [x*y*z]))
  @test ideal(Q, [2*x]) + ideal(Q, [x*(y+z)]) == ideal(Q, [x])
  @test iszero(ideal(Q, [y*z])*ideal(Q, [x]))

  a = quotient(ideal(Q, [zero(Q)]), ideal(Q, [y*z]))
  @test a == ideal(Q, [x])
  @test a == ideal(Q, gens(a))

  b = a + ideal(Q, [z])
  @test b == ideal(Q, [z, x])
  @test b == ideal(Q, gens(b))

  b = a*ideal(Q, [z])
  @test b == ideal(Q, [z*x])
  @test b == ideal(Q, gens(b))
end
