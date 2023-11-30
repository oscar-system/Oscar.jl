@testset "strands of graded modules" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])

  A = S[x y z; y z x]

  I = ideal(S, minors(A, 2))

  S1 = graded_free_module(S, 1)
  M, _ = quo(S1, (I*S1)[1])
  res = Oscar.SimpleComplexWrapper(free_resolution(M).C)[0:3]

  res_2 = Oscar.StrandComplex(res, 2)
  @test rank(res_2[0]) == 6
  @test rank(res_2[1]) == 3
  @test rank(res_2[2]) == 0
  @test !iszero(map(res_2, 1, (1,)))

  res_5, inc = Oscar.strand(res, 5)
  @test rank(res_5[0]) == 21
  @test rank(res_5[1]) == 30
  @test rank(res_5[2]) == 12
  @test !iszero(map(res_5, 1, (1,)))
  @test !iszero(map(res_5, 1, (2,)))
  @test iszero(compose(map(res_5, 1, (2,)), map(res_5, 1, (1,))))

  # Check commutativity of the squares
  for i in 0:1
      @test compose(map(res_5, 1, (i+1,)), inc[i]) == compose(inc[i+1], map(res, 1, (i+1,)))
  end
  @test Oscar.inclusion_map(res_5) === inc
end
