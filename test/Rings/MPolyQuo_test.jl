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

@testset "MpolyQuo.manipulation" begin
  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  I = ideal(R, [zero(R)])
  Q, q = quo(R,I)
  f = q(x*y)
  @test divides(one(Q), f) == (false, one(Q))

  A, _ = quo(R, 2*x^2-5*y^3)
  (x, y) = (A(x), A(y))

  @test iszero(x-x)
  @test x*deepcopy(x) == x^2
  @test iszero(det(matrix(A, [x 5; y^3 2*x])))

  @test !divides(x, y)[1]
  @test divides(x, x) == (true, one(A))
  @test divides(zero(A), x) == (true, zero(A))

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

  I = ideal(Q, [x^2*y-x+y,y+1])
  simplify!(I)
  @test I.SI[1] == singular_ring(Q)(-x+y) && I.SI[2] == singular_ring(Q)(y+1)
  J = ideal(Q, [x+y+1,y+1])
  @test issubset(J, I) == true
  @test issubset(I, J) == false
  @test (I == J) == false
  @test dim(J)  == 1
  @test dim(J)  == J.dim  # test case if dim(J) is already set
end
