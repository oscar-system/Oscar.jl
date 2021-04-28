@testset "mpoly_parametric.gcd" begin

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

@testset "mpoly_parametric.factor" begin

  function check_factor(a, esum)
    f = factor(a)

    @test isunit(unit(f))
    @test a == unit(f) * prod([p^e for (p, e) in f])
    @test esum == sum(e for (p, e) in f)

    f = factor_squarefree(a)

    @test isunit(unit(f))
    @test a == unit(f) * prod([p^e for (p, e) in f])
  end

  for (PR, Q) in ((PolynomialRing, QQ),
                  (PolynomialRing, GF(101)),
                  #(Singular.PolynomialRing, Singular.QQ) TODO enable after bug fixes
                 )

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
