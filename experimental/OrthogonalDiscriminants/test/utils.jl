@testset "order_omega_mod_N" begin
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 168) == (false, false)
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(6, 2, 168) == (true, false)
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 10) == (false, true)
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 12) == (true, true)
  for (d, q) in [(4, 2), (6, 3), (12, 5), (14, 7)]
    @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(d, q, 60)[2]
    @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(d, q, ZZ(60))[2]
  end
end

@testset "reduce_mod_squares" begin
  F, z = cyclotomic_field(7)
  for val in [F(0), F(1), -F(5), 1 + 2*z + z^5]
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(val) == val
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(12 * val) == 3 * val
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(27 * val) == 3 * val
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(val // 12) == 3 * val
  end
end

@testset "is_orthogonally_stable" begin
  t = character_table("J1");
  m = mod(t, 2);
  @test [is_orthogonally_stable(restrict(x, m)) for x in t] ==
        Bool[0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]
  @test [is_orthogonally_stable(restrict(x, m), check = false) for x in t] ==
        Bool[0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]
end

@testset "_orthogonal_discriminant_indicator0" begin
  t = character_table("L3(2)")
  @test Oscar.OrthogonalDiscriminants._orthogonal_discriminant_indicator0(t[2]) == "-7"
  @test Oscar.OrthogonalDiscriminants._orthogonal_discriminant_indicator0(mod(t, 3)[2]) == "O-"
  @test_throws ArgumentError Oscar.OrthogonalDiscriminants._orthogonal_discriminant_indicator0(t[4])
end

@testset "orbits under automorphisms" begin
  # _row_orbits_from_column_orbits
  M = identity_matrix(QQ, 5)
  M = M = [[M[i,j] for j in 1:5] for i in 1:5]
  orbs = [[1], [2, 3], [4, 5]]
  @test Oscar.OrthogonalDiscriminants._row_orbits_from_column_orbits(M, orbs) == orbs
  M = M[[2, 1, 4, 5, 3]]
  @test Oscar.OrthogonalDiscriminants._row_orbits_from_column_orbits(M, orbs) ==        [[1, 5], [2], [3, 4]]

  # _common_orbits
  @test Oscar.OrthogonalDiscriminants._common_orbits(
        [[1, 2], [3, 4]], [[1], [2, 3], [4]]) == [[1, 2, 3, 4]]
  @test Oscar.OrthogonalDiscriminants._common_orbits(
        [[1], [2, 3], [4]], [[1, 2], [3, 4]]) == [[1, 2, 3, 4]]
  @test Oscar.OrthogonalDiscriminants._common_orbits(
        [[1], [2], [3, 4], [5, 6]], [[1, 3], [2, 6], [4], [5]]) ==
        [[1, 3, 4], [2, 5, 6]]
  @test Oscar.OrthogonalDiscriminants._common_orbits(
        [[1, 3], [2, 6], [4], [5]], [[1], [2], [3, 4], [5, 6]]) ==
        [[1, 3, 4], [2, 5, 6]]

  # _inverse_fusion
  @test Oscar.OrthogonalDiscriminants._inverse_fusion(Int[], 1) == [[]]
  @test Oscar.OrthogonalDiscriminants._inverse_fusion([1, 2, 3, 3], 4) ==
        [[1], [2], [3, 4], []]
  @test Oscar.OrthogonalDiscriminants._inverse_fusion([1, 3, 3], 4) ==
        [[1], [], [2, 3], []]

  # orbits_group_automorphisms for ordinary tables
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("A5")) == [[1], [2, 3], [4], [5]]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("A6")) == [[1], [2, 3], [4, 5], [6], [7]]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("L3(8)")) == [[1], [2], [3, 4, 5, 6, 7, 8],
        [9, 10, 11, 12, 13, 14], [15, 16, 17, 18, 19, 20],
        [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, 32],
        [33], [34, 35, 36], [37, 38, 39, 40, 41, 42],
        [43, 44, 45, 46, 47, 48], [49, 50, 51, 52, 53, 54],
        [55, 56, 57, 58, 59, 60], [61], [62, 63, 64, 65, 66, 67],
        [68, 69], [70, 71, 72]]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("M11")) == [[x] for x in 1:10]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("U3(5)")) == [[1], [2], [3], [4, 5, 6], [7], [8],
        [9], [10], [11, 12], [13, 14]]

  # orbits_group_automorphisms for Brauer tables
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("A5", 2)) == [[1], [2, 3], [4]]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("A6", 2)) == [[1], [2, 3], [4, 5]]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("L3(8)", 2)) == [[1], [2, 3, 4, 5, 6, 7],
        [8, 9, 10], [11, 12, 13, 14, 15, 16], [17, 18, 19, 20, 21, 22],
        [23, 24, 25, 26, 27, 28], [29, 30, 31, 32, 33, 34],
        [35, 36], [37, 38, 39, 40, 41, 42], [43, 44, 45],
        [46, 47, 48, 49, 50, 51], [52, 53, 54, 55, 56, 57],
        [58, 59, 60, 61, 62, 63], [64]]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("M11", 2)) == [[x] for x in 1:5]
  @test Oscar.OrthogonalDiscriminants.orbits_group_automorphisms(
        character_table("U3(5)", 2)) == [[1], [2], [3, 4, 5], [6],
        [7, 8]]
end
