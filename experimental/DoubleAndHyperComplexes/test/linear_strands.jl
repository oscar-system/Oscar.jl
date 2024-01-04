@testset "linear strands of complexes of graded modules" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])

  S1 = graded_free_module(S, [0])
  I, inc = sub(S1, [f*S1[1] for f in gens(S)])
  M = cokernel(inc)

  res, aug = free_resolution(Oscar.SimpleFreeResolution, M)

  str, inc_str = Oscar.linear_strand(res)

  @test rank(str[0]) == 1
  @test rank(str[1]) == 3
  @test rank(str[2]) == 3
  @test rank(str[3]) == 1
  @test iszero(str[4])
  @test !Oscar.can_compute_index(str, (-1,))

  @test compose(inc_str[(1,)], map(res, 1, (1,))) == compose(map(str, 1, (1,)), inc_str[(0,)])
  @test compose(inc_str[(2,)], map(res, 1, (2,))) == compose(map(str, 1, (2,)), inc_str[(1,)])
  @test compose(inc_str[(3,)], map(res, 1, (3,))) == compose(map(str, 1, (3,)), inc_str[(2,)])

  K, pr = cokernel(inc_str)
  @test all(k->iszero(K[k]), 1:4)


  res_shift = shift(res, 1)

  str, inc_str = Oscar.linear_strand(res_shift)
  @test all(k->iszero(str[(k,)]), 0:10)

  S2 = graded_free_module(S, [0, 1])
  I1, inc1 = sub(S2, [f*S2[1] for f in gens(S)])
  I2, inc2 = sub(S2, [f*S2[2] for f in gens(S)])

  I = I1 + I2
  M, pr = quo(S2, I)

  res, aug = free_resolution(Oscar.SimpleFreeResolution, M)

  str, inc_str = Oscar.linear_strand(res)

  @test rank(str[0]) == 1
  @test rank(str[1]) == 3
  @test rank(str[2]) == 3
  @test rank(str[3]) == 1
  @test iszero(str[4])
  @test !Oscar.can_compute_index(str, (-1,))

  @test compose(inc_str[(1,)], map(res, 1, (1,))) == compose(map(str, 1, (1,)), inc_str[(0,)])
  @test compose(inc_str[(2,)], map(res, 1, (2,))) == compose(map(str, 1, (2,)), inc_str[(1,)])
  @test compose(inc_str[(3,)], map(res, 1, (3,))) == compose(map(str, 1, (3,)), inc_str[(2,)])

  K, pr = cokernel(inc_str)
  @test rank(K[0]) == 1
  @test rank(K[1]) == 3
  @test rank(K[2]) == 3
  @test rank(K[3]) == 1
  @test iszero(K[4])

  @test all(k->iszero(compose(inc_str[(k,)], pr[(k,)])), 0:3)
  @test all(k->compose(pr[(k,)], map(K, 1, (k,))) == compose(map(res, 1, (k,)), pr[(k-1,)]), 1:3)
  
  # Now for classical free resolutions
  res = free_resolution(M)

  str, inc_str = Oscar.linear_strand(res)

  @test rank(str[0]) == 1
  @test rank(str[1]) == 3
  @test rank(str[2]) == 3
  @test rank(str[3]) == 1
  @test iszero(str[4])
  @test !Oscar.can_compute_index(str, (-1,))

  @test compose(inc_str[(1,)], map(res, 1)) == compose(map(str, 1), inc_str[(0,)])
  @test compose(inc_str[(2,)], map(res, 2)) == compose(map(str, 2), inc_str[(1,)])
  @test compose(inc_str[(3,)], map(res, 3)) == compose(map(str, 3), inc_str[(2,)])

  K, pr = cokernel(inc_str)
  @test rank(K[0]) == 1
  @test rank(K[1]) == 3
  @test rank(K[2]) == 3
  @test rank(K[3]) == 1
  @test iszero(K[4])

  @test all(k->iszero(compose(inc_str[(k,)], pr[(k,)])), 0:3)
  @test all(k->compose(pr[(k,)], map(K, 1, (k,))) == compose(map(res, k), pr[(k-1,)]), 1:3)
  
  
  # Test linear strands of ComplexOfMorphism, too.
  res = free_resolution(M)

  str, inc_str = Oscar.linear_strand(res.C[3:-1:0])
  @test_throws ErrorException str[-1] # Not a free module in degree -1
end
