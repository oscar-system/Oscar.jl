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

end
