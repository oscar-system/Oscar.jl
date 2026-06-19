@testset "Basic modular group operations" begin
  s = @perm (2, 3)
  t = @perm (1, 2)
  S = symmetric_group(3)
  G = modular_subgroup_via_right_action(S(s), S(t))

  @test (@inferred index(G)) == 3

  s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
  t = cperm([1,4], [2,5,9,10,8], [3,7,6])
  G = modular_subgroup_via_right_action(s, t)
  @test (@inferred r_right_action(G)) == cperm([1,7,9,10,6], [2,3], [4,5,8])
  @test (@inferred j_right_action(G)) == cperm([1,8,3], [2,4,6], [5,7,10])

  SL2Z, S, T = Oscar._SL2Z_fp()
  actual_gens = @inferred word_gens(G)
  expected_gens = [
    S^-2,
    T^-2,
    S*T*S*T^-1*S^-1*T*S^-1,
    S*T^2*S*T^2*S^-1,
    S*T^5*S^-1
  ]
  @test Set(actual_gens) == Set(expected_gens)

  actual_mat_gens = gens(G)
  expected_mat_gens = [[-1 0; 0 -1],
    [1 -2; 0 1],
    [2 -1; -3 2],
    [2 -1; -3 2],
    [1 0; -5 1]
  ]
  @test Set(actual_mat_gens) == Set([ZZMatrix(m) for m in expected_mat_gens])

  M = matrix(ZZ, [1 0; -4 1])
  @test (@inferred s_t_decomposition(M)) == S * T^4 * S^-1
end
