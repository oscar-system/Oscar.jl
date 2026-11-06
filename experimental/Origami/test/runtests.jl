@testset "Basic origami operations" begin
  h = @perm (1, 2)
  v = @perm (2, 3)
  o = origami(h, v)

  @test degree(o) == max(degree(h), degree(v))
  @test (@inferred degree(o)) == 3

  G = perm_group(o)

  @test horizontal_perm(o) == G(h)
  @test vertical_perm(o) == G(v)
  @test stratum(o) == [2]
  @test genus(o) == 2
  @test index_monodromy_group(o) == 1
  @test cylinder_structure(o) == [(1, 2), (1, 1)]

  h2 = @perm (1, 2)(3, 4)
  v2 = @perm (2, 5, 3, 6)(4, 7)
  o2 = origami(h2, v2)

  @test (@inferred stratum(o2)) == [2, 1, 1]
  @test (@inferred genus(o2)) == 3
  @test (@inferred index_monodromy_group(o2)) == 2
  @test (@inferred sum_of_lyapunov_exponents(o2)) == 11//6
  @test (@inferred cylinder_structure(o2)) == [(1, 2), (1, 2), (1, 1), (1, 1), (1, 1)]

  h3 = @perm (1, 2, 3, 4, 5, 6)
  v3 = @perm (1, 6, 5, 4, 3, 2)
  o3 = origami(h3, v3)

  @test issetequal(
    translations(o3),
    [one(perm_group(o3)), cperm([1, 2, 3, 4, 5, 6]),
      cperm([1, 3, 5], [2, 4, 6]), cperm([1, 4], [2, 5], [3, 6]),
      cperm([1, 5, 3], [2, 6, 4]),
      cperm([1, 6, 5, 4, 3, 2])],
  )
  @test (@inferred is_hyperelliptic(o3)) == true

  h4 = @perm (1, 2, 3, 4, 5, 6)
  v4 = @perm (1, 4)(2, 6, 5, 3)
  o4 = origami(h4, v4)
  @test (@inferred translations(o4)) ==
        [one(perm_group(o4)), cperm([1, 4], [2, 5], [3, 6])] ||
    (@inferred translations(
    o4
  )) == [
    cperm([1, 4], [2, 5], [3, 6]), one(perm_group(o4))
  ]
  @test (@inferred is_hyperelliptic(o4)) == false
  @test (@inferred point_reflections(o4)) ==
    [cperm([2, 6], [3, 5]), cperm([1, 4], [2, 3], [5, 6])]

  @test (@inferred are_equivalent(o3, o4)) == false
  v5 = @perm (2, 5)(1, 6, 4, 3)
  o5 = origami(h4, v5)
  @test are_equivalent(o4, o5) == true
end

@testset "SL(2,Z) action" begin
  o = origami(cperm([1, 6, 4, 7, 5, 3], [2, 8]), cperm([1, 4, 5, 3, 8, 2, 6]))
  S = matrix(ZZ, [0 -1; 1 0])
  T = matrix(ZZ, [1 1; 0 1])

  @test action_sl2(identity_matrix(ZZ, 2), o) == o
  @test action_sl2(S, o) == action_s(o)
  @test action_sl2(T, o) == action_t(o)
  @test action_sl2(inv(S), o) == action_s_inv(o)
  @test action_sl2(inv(T), o) == action_t_inv(o)
  @test action_sl2(-one(S), o) == action_s(action_s(o))
  # The doctest matrix is the word S*T, applied from the right.
  @test action_sl2(ZZ[0 -1; 1 1], o) == action_s(action_t(o))
  @test action_sl2(ZZ[0 -1; 1 1], o) ==
        origami(cperm([1, 4, 6, 8], [5, 7]), cperm([1, 6, 4, 7, 5, 3], [2, 8]))

  # Agree with GAP's ActionOfSL2, including quotients that truncate toward zero.
  for A in (
      ZZ[0 -1; 1 1],
      matrix(ZZ, [1 0; -4 1]),
      matrix(ZZ, [2 1; 1 1]),
      matrix(ZZ, [-3 2; -2 1]),
      matrix(ZZ, [5 -3; 2 -1]),
      matrix(ZZ, [1 0; -3 1]),
      S^3 * T^-5 * S * T^2,
    )
    gap_o = GAP.Globals.ActionOfSL2(GapObj(A), GapObj(o))
    d = GAP.Globals.DegreeOrigami(gap_o)::Int
    G = symmetric_group(d)
    h = PermGroupElem(G, GAP.Globals.HorizontalPerm(gap_o))
    v = PermGroupElem(G, GAP.Globals.VerticalPerm(gap_o))
    @test action_sl2(A, o) == origami(h, v; check=false)
  end

  # T^N for an exponent that does not fit in an Int. The decomposition of
  # [1 0; -N 1] is the word S*T^N*S^{-1}.
  N = ZZ(2)^100
  o1 = action_s_inv(o)
  h = horizontal_perm(o1)
  steps = Int(mod(N, order(Int, h)))
  expected = o1
  for _ in 1:steps
    expected = action_t(expected)
  end
  expected = action_s(expected)
  @test action_sl2(matrix(ZZ, [1 0; -N 1]), o) == expected

  @test_throws ArgumentError action_sl2(matrix(ZZ, [1 0; 0 -1]), o)
  @test_throws ArgumentError action_sl2(identity_matrix(ZZ, 3), o)
end

@testset "Listing strata" begin
  include("stratum_lists.jl")
  actual_strat_list = origamis([1, 1, 1, 1], 8)
  actual_strat_list_normal = [normal_form(o) for o in actual_strat_list]
  expected_strat_list_normal = [normal_form(o) for o in expected_strat_list]
  @test issetequal(expected_strat_list_normal, actual_strat_list_normal)
end

@testset "Action of S on the extended homology of T_n" begin
  # Check the S action on the extended homology of T_3 against its known matrix and
  # against the action on the non-extended homology, which it contains as the upper left
  # block, followed by the image of the additional loop and a zero row
  M = action_of_s_on_extended_homology_of_tn(3)
  expected = matrix(ZZ, [
     0 -1  0  0  0  0  0  0  0  0  0;
     1  0  0  0  0  0  0  0  0  0  0;
    -1  0  1  0  0  0  0 -1  0  0  0;
     0  0  0  0  0  0  0 -1  1  0  0;
     0  0  0  0  0  1  0 -1  0  0  0;
    -1  0  0  1  0  0  0 -1  0  0  0;
     0  0  0  0  0  0  0 -1  0  1  0;
     0  0  0  0  0  0  1 -1  0  0  0;
    -1  0  0  0  1  0  0 -1  0  0  0;
     0  0  0  0  0  0  0 -1  0  0  1;
     0  0  0  0  0  0  0  0  0  0  0
  ])
  @test M == expected

  A = action_of_s_on_homology_of_tn(3)
  N = nrows(A)
  c = matrix(ZZ, N, 1, [0; 0; fill(-1, N - 2)])
  @test M[1:N, 1:N] == A
  @test M[1:N, N + 1:N + 1] == A * c
  @test iszero(M[N + 1:N + 1, :])
end

@testset "Action of T on the extended homology of T_n" begin
  # Check the T action on the extended homology of T_3 against its known matrix and
  # against the action on the non-extended homology, which it contains as the upper left
  # block, followed by the image of the additional loop and a zero row
  M = action_of_t_on_extended_homology_of_tn(3)
  expected = matrix(ZZ, [
    1  1  0  0  0  0  0  0  0  0  0;
    0  1  0  0  0  0  0  0  0  0  0;
    0  0  1  0  0  0  0  0 -1  0  0;
    0  0  0  1  0  0  0  0 -1  0  0;
    0  0  0  0  1  0  0  0 -1  0  0;
    0 -1  0  0  0  0  0  1 -1  0  0;
    0  0  0  0  0  1  0  0 -1  0  0;
    0 -1  0  0  0  0  1  0 -1  0  0;
    0 -1  0  0  0  0  0  0 -1  1  0;
    0  0  0  0  0  0  0  0 -1  0  1;
    0  0  0  0  0  0  0  0  0  0  0
  ])
  @test M == expected

  A = action_of_t_on_homology_of_tn(3)
  N = nrows(A)
  c = matrix(ZZ, N, 1, [0; 0; fill(-1, N - 2)])
  @test M[1:N, 1:N] == A
  @test M[1:N, N + 1:N + 1] == A * c
  @test iszero(M[N + 1:N + 1, :])
end

@testset "Relations of SL_2(Z) acting on the homology of T_n" begin
  # Check that the actions of S and T on the homology of T_n satisfy the relations of
  # SL_2(Z), and that S has order 4
  for n in 2:3
    S = action_of_s_on_homology_of_tn(n)
    T = action_of_t_on_homology_of_tn(n)
    I = identity_matrix(ZZ, n^2 + 1)
    @test S^4 == I
    @test S^2 != I
    @test (S * T)^3 == S^2
  end
end

@testset "Relations of SL_2(Z) acting on the extended homology of T_n" begin
  # Check that the actions on the extended homology of T_n satisfy the relations of SL_2(Z)
  # between non-empty words, that S^4 is the projection EF, that the sum of all puncture loops
  # is mapped to zero, and that action_of_matrix_on_extended_homology_of_tn agrees with S and T
  for n in 2:3
    S = action_of_s_on_extended_homology_of_tn(n)
    T = action_of_t_on_extended_homology_of_tn(n)
    @test (S * T)^3 == S^2
    @test S^5 == S
    @test S^8 == S^4
    @test S^4 != identity_matrix(ZZ, n^2 + 2)
    N = n^2 + 1
    c = matrix(ZZ, N, 1, [0; 0; fill(-1, N - 2)])
    @test S^4 == vcat(hcat(identity_matrix(ZZ, N), c), zero_matrix(ZZ, 1, N + 1))
    r = matrix(ZZ, n^2 + 2, 1, [0; 0; fill(1, n^2)])
    @test iszero(S * r)
    @test iszero(T * r)
    @test action_of_matrix_on_extended_homology_of_tn(n, [0 -1; 1 0]) == S
    @test action_of_matrix_on_extended_homology_of_tn(n, [1 1; 0 1]) == T
    @test action_of_matrix_on_extended_homology_of_tn(n, [0 -1; 1 1]) == S * T
  end
end
