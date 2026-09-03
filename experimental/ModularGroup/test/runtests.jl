# Construction of modular subgroups from coset actions, their generators,
# membership tests and the S/T decomposition of SL(2,Z) matrices.
@testset "Basic modular group operations" begin
  s = @perm (2, 3)
  t = @perm (1, 2)
  S = symmetric_group(3)
  G = modular_subgroup_via_right_action(S(s), S(t))

  @test (@inferred index(G)) == 3

  s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
  t = cperm([1,4], [2,5,9,10,8], [3,7,6])
  G = modular_subgroup_via_right_action(s, t)
  @test (@inferred r_right_perm(G)) == cperm([1,7,9,10,6], [2,3], [4,5,8])
  @test (@inferred j_right_perm(G)) == cperm([1,8,3], [2,4,6], [5,7,10])

  SL2Z = Oscar._SL2Z_fp()
  S, T = gens(SL2Z)
  actual_gens = @inferred word_gens(G)
  expected_gens = [
    S^-2,
    T^-2,
    S*T*S*T^-1*S^-1*T*S^-1,
    S*T^2*S*T^2*S^-1,
    S*T^5*S^-1
  ]
  @test issetequal(actual_gens, expected_gens)

  actual_mat_gens = @inferred gens(G)
  expected_mat_gens = [[-1 0; 0 -1],
    [1 -2; 0 1],
    [2 -1; -3 2],
    [1 0; -5 1]
  ]
  @test issetequal(matrix.(actual_mat_gens), [ZZMatrix(m) for m in expected_mat_gens])
  @test ZZMatrix(expected_mat_gens[3]) in G
  @test all(g in G for g in actual_mat_gens)

  M = matrix(ZZ, [1 0; -4 1])
  @test (@inferred s_t_decomposition(M)) == S * T^4 * S^-1
  w = expected_gens[3]
  @test (@inferred Oscar._image_of_pt(w, G, 1)) == 1
  @test is_word_elm_of(w, G)
end

@testset "Modular group equality, hashing and membership" begin
  s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
  t = cperm([1,4], [2,5,9,10,8], [3,7,6])
  G = modular_subgroup_via_right_action(s, t)

  # relabelling the cosets (fixing coset 1) yields the same subgroup
  c = cperm([2,3], [5,8,10])
  G2 = modular_subgroup_via_right_action(s^c, t^c)
  @test G == G2
  @test issubset(G, G2)
  @test hash(G) == hash(G2)
  @test length(Set([G, G2])) == 1

  # elements with equal matrices are equal and hash equally
  A = matrix(ZZ, [1 -2; 0 1])
  x = Oscar.ModularGroupElem(G, A)
  y = Oscar.ModularGroupElem(G, deepcopy(A))
  @test x == y
  @test hash(x) == hash(y)
  @test length(Set([x, y])) == 1
  @test x in G
  @test one(G) in G

  # permutations from different symmetric groups are coerced
  H = modular_subgroup_via_right_action(cperm([2,3]), cperm([1,2]))
  @test index(H) == 3
  @test H != G
  S3 = symmetric_group(3)
  @test H == modular_subgroup_via_right_action(S3(cperm([2,3])), S3(cperm([1,2])))

  # the whole group SL(2,Z)
  W = modular_subgroup_via_right_action(cperm(), cperm())
  @test index(W) == 1
  @test matrix(ZZ, [1 0; 1 1]) in W

  # only 2x2 matrices of determinant 1 can be elements
  @test !(identity_matrix(ZZ, 3) in G)
  @test !(matrix(ZZ, [1 0; 0 -1]) in G)
  @test_throws ArgumentError s_t_decomposition(identity_matrix(ZZ, 3))
  @test_throws ArgumentError s_t_decomposition(matrix(ZZ, [1 0; 0 -1]))

  # membership must not be linear in the size of the matrix entries;
  # T^N lies in G iff N is a multiple of the length of the t-cycle of 1
  N = ZZ(10)^12
  @test matrix(ZZ, [1 N; 0 1]) in G
  @test !(matrix(ZZ, [1 N + 1; 0 1]) in G)
end
