@testset "mpoly_affine_algebras.normalization" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  Q, _ = quo(R, ideal(R, [z - x^4, z - y^6]))
  for (S, M, FI) in normalization(Q)
    @test parent(FI[1]) == Q
    @test isa(FI[2], Oscar.Ideal)
    @test parent(M(Q(x+y))) == S
  end

  stuff, deltas, tot_delta = normalization_with_delta(Q)
  @test isa(tot_delta, Int)
  @test length(stuff) == length(deltas)
  for ((S, M, FI), delta) in zip(stuff, deltas)
    @test parent(FI[1]) == Q
    @test isa(FI[2], Oscar.Ideal)
    @test parent(M(Q(x+y))) == S
    @test isa(delta, Int)
  end

  for algorithm in (:equidimDec, :primeDec)
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
    Q, _ = quo(R, ideal(R, [(x^2 - y^3)*(x^2 + y^2)*x]))
    for (S, M, FI) in normalization(Q, alg=algorithm)
      @test parent(FI[1]) == Q
      @test isa(FI[2], Oscar.Ideal)
      @test parent(M(Q(x+y))) == S
    end

    stuff, deltas, tot_delta = normalization_with_delta(Q, alg=algorithm)
    @test isa(tot_delta, Int)
    @test length(stuff) == length(deltas)
    for ((S, M, FI), delta) in zip(stuff, deltas)
      @test parent(FI[1]) == Q
      @test isa(FI[2], Oscar.Ideal)
      @test parent(M(Q(x+y))) == S
      @test isa(delta, Int)
    end
  end

  @test !isnormal(quo(R, ideal(R, [x^2 - y^3]))[1])
  @test isnormal(quo(R, ideal(R, [x - y^3]))[1])

  R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
  @test_throws ArgumentError isnormal(quo(R, ideal(R, [x - y^3]))[1])
end

@testset "mpoly_affine_algebras.integral_basis" begin

  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2)
  @test all(p -> base_ring(parent(p)) == R, nums)
  @test base_ring(parent(den)) == R
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 3)
  @test_throws ArgumentError integral_basis((x+y)*(x+y^2), 1)

  R, (x, y) = Singular.PolynomialRing(Singular.QQ, ["x", "y"])
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2)
  @test all(p -> base_ring(parent(p)) == R, nums)
  @test base_ring(parent(den)) == R
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 3)
  @test_throws ArgumentError integral_basis((x+y)*(x+y^2), 1)

  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 2)

  R, (x, y) = PolynomialRing(GF(2), ["x", "y"])
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)
end

@testset "Noether normalization" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  A, _ = quo(R, ideal(R, [x*y, x*z]))
  L = noether_normalization(A)
  @test length(L) == 3 # just a smoke test
end
