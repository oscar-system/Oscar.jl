begin

local check_gcd, check_factor

function check_gcd(a, b, gdiv)
  g = gcd(a, b)
  if iszero(g)
    @test iszero(a)
    @test iszero(b)
    return
  end
  @test iszero(gdiv) || divides(g, gdiv)[1]
  a = divexact(a, g)
  b = divexact(b, g)
  @test isunit(gcd(a, b))
end

function check_factor(a, esum)
  f = factor(a)

  @test isunit(unit(f))
  @test a == unit(f) * prod([p^e for (p, e) in f])
  @test esum == sum(e for (p, e) in f)

  f = factor_squarefree(a)

  @test isunit(unit(f))
  @test a == unit(f) * prod([p^e for (p, e) in f])
end


@testset "mpoly-nested.frac.gcd" begin
  r, (t1, t2, t3) = PolynomialRing(QQ, ["t1", "t2", "t3"])
  r = FractionField(r)
  (t1, t2, t3) = map(r, (t1, t2, t3))
  r, (x1, x2, x3, x4) = PolynomialRing(r, ["x1", "x2", "x3", "x4"])

  check_gcd(zero(r), zero(r), zero(r))
  check_gcd(zero(r), x1+t1, x1+t1)
  check_gcd(x2+t2, zero(r), x2+t2)
  check_gcd(zero(r), zero(r), zero(r))

  g = (t2*x1+t1//(t2+t3)*x2+t3+x4*t1*t2*t3)*x1
  a = (x3*t2//t1+t3)*x2
  b = (x1+x3*x2*t1//(t2*t3))*x3
  check_gcd(g*a, g*b, g)
end

@testset "mpoly-nested.Q(y)[x]" begin
  Qx, x = PolynomialRing(QQ, "x")
  Fx = FractionField(Qx)
  x = Fx(x)
  R, y = PolynomialRing(Fx, "y")
  check_factor((y*inv(x)+x)^2, 2)
  @test isunit(gcd(x+y, x-y))

  Qx, (x,) = PolynomialRing(QQ, ["x"])
  Fx = FractionField(Qx)
  x = Fx(x)
  R, y = PolynomialRing(Fx, "y")
  check_factor((y*inv(x)+x)^2, 2)
  @test isunit(gcd(x+y, x-y))

  if false  # enable only for long tests
    Qx, x = PolynomialRing(QQ, "x")
    Fx = FractionField(Qx)
    x = Fx(x)
    R, (y, ) = PolynomialRing(Fx, ["y"])
    check_factor((y*inv(x)+x)^2, 2)
    @test isunit(gcd(x+y, x-y))

    Qx, (x,) = PolynomialRing(QQ, ["x"])
    Fx = FractionField(Qx)
    x = Fx(x)
    R, (y, ) = PolynomialRing(Fx, ["y"])
    check_factor((y*inv(x)+x)^2, 2)
    @test isunit(gcd(x+y, x-y))
  end
end


@testset "mpoly-nested.frac.factor" begin

  if false  # enable only for long tests
    rings = ((PolynomialRing, QQ),
             (PolynomialRing, GF(101)),
             (Singular.PolynomialRing, Singular.QQ)
            )
  else
    rings = ((PolynomialRing, GF(101)),)
  end

  for (PR, Q) in rings
    r, (t1, t2, t3) = PR(Q, ["t1", "t2", "t3"])
    r = FractionField(r)
    (t1, t2, t3) = map(r, (t1, t2, t3))
    r, (x1, x2, x3, x4) = PolynomialRing(r, ["x1", "x2", "x3", "x4"])

    check_factor((t1+t2*t3)*(x1*1//t2^2+x2*1//t2+x3*1//t3)*
                 (x1^2*t1+x2^2*t1//t2+x3^2*t1//t3), 2)

    check_factor(1//(t1+t2*t3)*(x1*1//t2^2+x2*1//t2+x3*1//t3+x4)^3*
                 (x1^2*t1+x2^2*t1//t2+x4*x3^2*t1//t3)^2, 5)
  end
end

@testset "mpoly-nested.iterated.conversion" begin
  R = QQ
  R, x1 = PolynomialRing(R, "x1")
  p = (x1 + 2)^2
  @test p == renest(R, denest(denest(R), p))

  R = QQ
  R, (x1, x2) = PolynomialRing(R, ["x1", "x2"])
  p = (x1*x2 + x1 + x2^2)^2
  @test p == renest(R, denest(denest(R), p))

  if false  # enable only for long tests
    R = QQ
    R, x3 = PolynomialRing(R, "x3")
    R, x2 = PolynomialRing(R, "x2")
    R, x1 = PolynomialRing(R, "x1")
    p = (x1*x2*x3 + x1 + x2^2 + x3^3)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, (x3, x4) = PolynomialRing(R, ["x3", "x4"])
    R, x2 = PolynomialRing(R, "x2")
    R, x1 = PolynomialRing(R, "x1")
    p = (x1*x2*x3*x4 + x1 + x2^2 + x3^3 + x4^4)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, x4 = PolynomialRing(R, "x4")
    R, (x2, x3) = PolynomialRing(R, ["x2", "x3"])
    R, x1 = PolynomialRing(R, "x1")
    p = (x1*x2*x3*x4 + x1 + x2^2 + x3^3 + x4^4)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, x4 = PolynomialRing(R, "x4")
    R, x3 = PolynomialRing(R, "x3")
    R, (x1, x2) = PolynomialRing(R, ["x1", "x2"])
    p = (x1*x2*x3*x4 + x1 + x2^2 + x3^3 + x4^4)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, x5 = PolynomialRing(R, "x5")
    R, (x3, x4) = PolynomialRing(R, ["x3", "x4"])
    R, (x1, x2) = PolynomialRing(R, ["x1", "x2"])
    p = (x1*x2*x3*x4*x5 + x1 + x2^2 + x3^3 + x4^4 + x5^5)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, (x4, x5) = PolynomialRing(R, ["x5", "x6"])
    R, x3 = PolynomialRing(R, "x4")
    R, (x1, x2) = PolynomialRing(R, ["x1", "x2"])
    p = (x1*x2*x3*x4*x5 + x1 + x2^2 + x3^3 + x4^4 + x5^5)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, (x4, x5) = PolynomialRing(R, ["x4", "x5"])
    R, (x2, x3) = PolynomialRing(R, ["x2", "x3"])
    R, x1 = PolynomialRing(R, "x1")
    p = (x1*x2*x3*x4*x5 + x1 + x2^2 + x3^3 + x4^4 + x5^5)^2
    @test p == renest(R, denest(denest(R), p))

    R = QQ
    R, (x5, x6) = PolynomialRing(R, ["x5", "x6"])
    R, (x3, x4) = PolynomialRing(R, ["x3", "x4"])
    R, (x1, x2) = PolynomialRing(R, ["x1", "x2"])
    p = (x1*x2*x3*x4*x5*x6 + x1 + x2^2 + x3^3 + x4^4 + x5^5 + x6^6)^2
    @test p == renest(R, denest(denest(R), p))
  end
end

@testset "mpoly-nested.iterated.gcd_factor" begin
  r3 = PolynomialRing(QQ, "x3")[1]
  r2 = PolynomialRing(r3, "x3")[1]

  R1 = PolynomialRing(r2, "x3")[1]
  xs1 = [gen(R1), R1(gen(r2)), R1(gen(r3))]

  R2 = PolynomialRing(r3, ["x1", "x2"])[1]
  xs2 = [gen(R2, 1), gen(R2, 2), R2(gen(r3))]

  for ((x1, x2, x3), R) in ((xs1, R1), (xs2, R2))
    check_gcd(zero(R), zero(R), zero(R))
    check_gcd(zero(R), x2+x3, x2+x3)
    check_gcd(x1+x2, zero(R), x1+x2)
    check_gcd(zero(R), zero(R), zero(R))

    g = (x1 + 2*x2 + 3*x3)
    a = (x1^2 + 4*x2^2 + 5*x3^2)
    b = (x1^3 + 6*x2^3 + 7*x3^3)
    check_gcd(g*a, g*b, g)

    check_factor(g^2*a*b, 4)
  end
end

end
