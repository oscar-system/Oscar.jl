GAP.Packages.load("StandardFF")

@testset "standard_finite_field unit tests" begin
  function MyPoly(p, n)
    return defining_polynomial(standard_finite_field(p, n))
  end

  function GAPPoly(p, n)
    poly = GAP.Globals.StandardFiniteField(Int(p), n)
    poly = GAP.Globals.DefiningPolynomial(poly)
    poly = GAP.Globals.CoefficientsOfUnivariatePolynomial(poly)
    k = Nemo.Native.GF(p)
    poly = polynomial(k, map(k, Vector{GAP.FFE}(poly)))
  end

  function compare_poly(p, n)
    F = Nemo.Native.GF(p)
    F.(collect(coefficients(MyPoly(p, n)))) == F.(collect(coefficients(GAPPoly(p, n))))
  end

  small_test_primes = [3, 5, 7]
  # medium_test_primes = ZZ.([3433, 4073])


  @testset "characteristic 2 test GF(2^(2^$k))" for k in 3:5
    @test compare_poly(2, 2^k)
    @test compare_poly(ZZ(2), 2^k)
  end
  @testset "small primes GF($p^(2^$k))" for p in small_test_primes, k in 3:5
    @test compare_poly(p, 2^k)
    @test compare_poly(ZZ(p), 2^k)
  end
  @testset "small primes GF($p^($p^3))" for p in small_test_primes
    @test compare_poly(p, p^3)
    @test compare_poly(ZZ(p), p^3)
  end
  @testset "composite tests" begin
    @test compare_poly(2, (2^3) * (3^3))
    @test compare_poly(3, (2^3) * 3 * 5)
    @test compare_poly(5, 2 * (3^3) * 5)
    @test compare_poly(ZZ(2), (2^3) * (3^3))
    @test compare_poly(ZZ(3), (2^3) * 3 * 5)
    @test compare_poly(ZZ(5), 2 * (3^3) * 5)
  end
end
