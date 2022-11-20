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

  # promote rule
  K = GF(2)
  Kx, (x, y) = K["x", "y"]
  I = ideal(Kx, elem_type(Kx)[])
  Kx, = quo(Kx, I)
  x = Kx(x)
  @test K(2) * x == Kx(2) * x
end

@testset "MPolyQuo.ideals" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  Q, _ = quo(R, ideal(R, [x*y, x*z]))
  (x, y, z) = map(Q, (x, y, z))

  I = ideal(Q, [x^2, y^2, z^2])
  @test I isa Oscar.Ideal
  @test y^2 in I
  @test x*(x-y) in I
  @test !(x in I)
  @test !(y*z in I)

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
  SQ = singular_poly_ring(Q)
  @test I.SI[1] == SQ(-x+y) && I.SI[2] == SQ(y+1)
  J = ideal(Q, [x+y+1,y+1])
  @test issubset(J, I) == true
  @test issubset(I, J) == false
  @test (I == J) == false
  @test dim(J)  == 1
  @test dim(J)  == J.dim  # test case if dim(J) is already set

  R, (x, y) = grade(PolynomialRing(QQ, [ "x", "y"])[1], [ 1, 2 ])
  I = ideal(R, [ x*y ])
  Q, RtoQ = quo(R, I)
  J = ideal(Q, [ x^3 + x*y, y, x^2 + y ])
  @test minimal_generating_set(J) == [ Q(y), Q(x^2 + y) ]
  @test minimal_generating_set(ideal(Q, [ Q() ])) == elem_type(Q)[]

  # 1530
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  J = ideal(Q, [x^7, y^2])
  @test !(one(Q) in J)

  # allowing empty set of generators 
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, x*y)
  Q, _ = quo(R, I)
  J = ideal(Q, elem_type(Q)[])
  @test Q(x*y) in J
  @test !(Q(x) in J)
  J2 = ideal(Q, elem_type(R)[])
  @test Q(x*y) in J2
  @test !(Q(x) in J2)
end

@testset "prime ideals" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x])
  A, _ = quo(R, I)
  J = ideal(A, [A(y)])
  @test is_prime(J)
  J2 = ideal(A, [A(z*y)])
  @test !is_prime(J2)
end

@testset "modulus - MPAnyQuoRing MPAnyNonQuoRing" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x+y+z])
  A, _ = quo(R,I)
  @test modulus(R) == ideal(R,[zero(R)])
  @test modulus(A) == I
  U=MPolyComplementOfKPointIdeal(R,[0,0,0])
  Rl,_ = Localization(R,U)
  Il = Rl(I)
  Al, _ = quo(Rl, Il)
  @test modulus(Rl) == ideal(Rl,[zero(Rl)])
  @test modulus(Al) == Il
  U2=MPolyComplementOfPrimeIdeal(ideal(R,[x^2+1,y^2+1,z]))
  Rl2,_ = Localization(R,U2)
  Il2 = Rl2(I)
  Al2,_ = quo(Rl2,Il2)
  @test modulus(Rl2) == ideal(Rl2,[zero(Rl2)])
  @test modulus(Al2) == Il2
  U3=MPolyPowersOfElement(x+y)
  Rl3,_ = Localization(R,U3)
  Il3 = Rl3(I)
  Al3,_ = quo(Rl3,Il3)
  @test modulus(Rl3) == ideal(Rl3,[zero(Rl3)])
  @test modulus(Al3) == Il3
end

@testset "saturated ideal compatibility" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  I = ideal(R, [x+y+z])
  A, _ = quo(R,I)
  J = ideal(A, [x, y])
  @test z in saturated_ideal(J)
end

