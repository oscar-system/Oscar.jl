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
end
