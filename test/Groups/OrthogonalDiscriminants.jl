@testset "all_od_infos" begin
  iob = IOBuffer()
  show_OD_info("G2(3)", iob)
  str1 = String(take!(iob))

  all_entries = all_od_infos();
  @test all_entries isa Vector

  show_OD_info("G2(3)", iob)
  str2 = String(take!(iob))
  @test str1 == str2

  @test length(all_entries) == length(all_od_infos(is_simple)) +
                               length(all_od_infos(! is_simple))
  @test length(all_entries) == length(all_od_infos(is_sporadic_simple)) +
                               length(all_od_infos(! is_sporadic_simple))
  @test length(all_od_infos(is_sporadic_simple, is_simple)) ==
                               length(all_od_infos(is_sporadic_simple))
  @test length(all_od_infos(is_sporadic_simple, ! is_simple)) == 0
  A5_entries = all_od_infos(identifier => "A5")
  @test length(A5_entries) == 3
  @test length(all_od_infos(identifier => ["A5", "A6"])) ==
        length(A5_entries) +
        length(all_od_infos(identifier => "A6"))
  A_entries = all_od_infos(identifier => x -> x[1] == 'A')
  @test A5_entries == A_entries[1:3]
  @test A5_entries == all_od_infos(is_simple, identifier => "A5")
  @test all_od_infos(identifier => ["A5", "A6"], characteristic => [3, 5]) ==
        all_od_infos(identifier => ["A5", "A6"], characteristic => isodd)
  @test all_od_infos(characteristic => [0, 2]) ==
        all_od_infos(characteristic => iseven)
  @test length(all_od_infos(dim => isodd)) == 0
  A_dims = collect(Set([parse(Int, filter(isdigit, x[:charname])) for x in A_entries]))
  @test all_od_infos(identifier => x -> x[1] == 'A', dim => A_dims) == A_entries
  Oplusminus = all_od_infos(orthogonal_discriminant => ["O+", "O-"]);
  @test length(all_od_infos(orthogonal_discriminant => "O+")) +
        length(all_od_infos(orthogonal_discriminant => "O-")) ==
        length(Oplusminus)
  @test length(Oplusminus) <= length(all_od_infos(characteristic => is_positive))
  deg_1 = all_od_infos(degree => 1);
  deg_other = all_od_infos(degree => (x -> x > 1));
  @test length(all_entries) == length(deg_1) + length(deg_other)
  @test deg_1 == all_od_infos(degree => [1]);
end

@testset "orthogonal_discriminant(s)" begin
  for name in ["S4", "A5", "A6", "J1"]
    t = character_table(name)
    vals = orthogonal_discriminants(t)
    @test vals == [orthogonal_discriminant(chi) for chi in t]
    @test length(vals) == length(t)
    @test vals[1] == ""
    for p in vcat([i for (i, e) in factor(order(t))], [101])
      modt = mod(t, p)
      vals = orthogonal_discriminants(modt)
      @test vals == [orthogonal_discriminant(chi) for chi in modt]
      @test length(vals) == length(modt)
      @test vals[1] == ""
    end
  end
end

@testset "orthogonal_discriminants for central extensions" begin
  name = "2.L3(4)"
  t = character_table(name)
  @test orthogonal_discriminants(t) ==
  ["","21","","","","","","","","105","","","1","1","1","1","-1","-1"]
end


# Make sure that the tests start without cached character tables.
# (Running the tests will store information that changes some test outputs,
# thus running the tests twice needs these calls.)
GAP.Globals.UnloadCharacterTableData()
empty!(Oscar.character_tables_by_id)

Oscar.@_AuxDocTest "show and print character tables", (fix = false),
raw"""
a table with irrational ODs

```jldoctest show_with_ODs.test
julia> using Oscar

julia> t = character_table("L2(11)");

julia> Oscar.show_with_ODs(t, stdout)
L2(11)

                2  2  2  1  .  .  1   .   .
                3  1  1  1  .  .  1   .   .
                5  1  .  .  1  1  .   .   .
               11  1  .  .  .  .  .   1   1
                                           
                  1a 2a 3a 5a 5b 6a 11a 11b
               2P 1a 1a 3a 5b 5a 3a 11b 11a
               3P 1a 2a 1a 5b 5a 2a 11a 11b
               5P 1a 2a 3a 1a 1a 6a 11a 11b
              11P 1a 2a 3a 5a 5b 6a  1a  1a
    d      OD   2                          
X_1 1           +  1  1  1  1  1  1   1   1
X_2 2           o  5  1 -1  .  .  1   B  /B
X_3 2           o  5  1 -1  .  .  1  /B   B
X_4 1     -11   + 10 -2  1  .  .  1  -1  -1
X_5 1     -11   + 10  2  1  .  . -1  -1  -1
X_6 1           + 11 -1 -1  1  1 -1   .   .
X_7 2 22-11b5   + 12  .  .  A A*  .   1   1
X_8 2 11b5+33   + 12  .  . A*  A  .   1   1

A = -z_5^3 - z_5^2 - 1
A* = z_5^3 + z_5^2
B = z_11^9 + z_11^5 + z_11^4 + z_11^3 + z_11
/B = -z_11^9 - z_11^5 - z_11^4 - z_11^3 - z_11 - 1
```

a table with unknown ODs (outside the scope of the database)

```jldoctest show_with_ODs.test
julia> t = character_table("D8");

julia> Oscar.show_with_ODs(t, stdout)
D8

          2  3  3  2  2  2
                          
            1a 2a 4a 2b 2c
         2P 1a 1a 2a 1a 1a
    d OD  2               
X_1 1     +  1  1  1  1  1
X_2 1     +  1  1  1 -1 -1
X_3 1     +  1  1 -1  1 -1
X_4 1     +  1  1 -1 -1  1
X_5 1  ?  +  2 -2  .  .  .
```
"""

@testset "compute via the Atlas of Group Representations" begin
#TODO: How can we test functions that need web access?
#  # Check equality of positive results in small cases.
#  l = all_od_infos(comment_matches => "AGR", dim => 1:8);
#  for entry in l
#    chi = Oscar.character_of_entry(entry)
#    comp = Oscar.od_from_atlas_group(chi)
#    @test comp == (true, entry[:valuestring])
#  end
#
#  # Check negative results.
#  t = character_table("L2(17)")
#  chi = t[10]
#  @test Oscar.od_from_atlas_group(chi) == (false, "")
#  t = character_table("L2(13)", 3)
#  chi = t[4]
#  @test Oscar.od_from_atlas_group(chi) == (false, "")
end

@testset "dimensions and Gram determinants of Specht modules" begin
  @test [dimension_specht_module(p) for p in partitions(5)] ==
        [1, 4, 5, 6, 5, 4, 1]
  @test [gram_determinant_specht_module(p) for p in partitions(5)] ==
        [[],
         [[5, 1]],
         [[2, 1], [3, 4]],
         [[2, 6], [5, 3]],
         [[2, 14], [3, 1]],
         [[2, 4], [3, 4], [5, 3]],
         [[2, 3], [3, 1], [5, 1]]]
end

@testset "group order" begin
  d = 8
  q = 3
  order_oplus = order(omega_group(1, d, q))
  order_ominus = order(omega_group(-1, d, q))
  f = Oscar.order_omega_mod_N
  @test f(d, q, order_oplus) == (true, false)
  @test f(d, q, order_ominus) == (false, true)
  @test f(d, q, order(omega_group(0, d-1, q))) == (true, true)
  @test f(d, q, q*order_oplus) == (false, false)
  @test f(d, q, (q-1)*order_oplus) == (false, false)

  @test_throws ArgumentError f(5, 2, 1)

  for entry in all_od_infos(comment_matches => "order")
    chi = Oscar.character_of_entry(entry)
    @test Oscar.od_from_order(chi) == (true, entry[:valuestring])
  end
  for entry in all_od_infos(identifier => "A8")
    if ! comment_matches(entry, "order")
      chi = Oscar.character_of_entry(entry)
      @test Oscar.od_from_order(chi) == (false, "")
    end
  end
end

@testset "eigenvalues" begin
  for entry in all_od_infos(comment_matches => "ev")
    chi = Oscar.character_of_entry(entry)
    @test Oscar.od_from_eigenvalues(chi) == (true, entry[:valuestring])
  end
  for entry in all_od_infos(identifier => "A8")
    if ! comment_matches(entry, "ev")
      chi = Oscar.character_of_entry(entry)
      @test Oscar.od_from_eigenvalues(chi) == (false, "")
    end
  end
end

@testset "Specht modules" begin
  for entry in all_od_infos(comment_matches => "specht")
    chi = Oscar.character_of_entry(entry)
    @test Oscar.od_for_specht_module(chi) == (true, entry[:valuestring])
  end
  for entry in all_od_infos(identifier => "A8")
    if ! comment_matches(entry, "specht")
      chi = Oscar.character_of_entry(entry)
      @test Oscar.od_for_specht_module(chi) == (false, "")
    end
  end
end

@testset "p-groups" begin
  rx = r"syl\((?<p>\d+)\)"
  info = all_od_infos(characteristic => 0,
                 comment_matches => (x -> any(startswith("syl"), x)))
  for entry in info[1:100]
    chi = Oscar.character_of_entry(entry)
    c = filter(startswith("syl"), entry[:comment])
    for e in c
      p = parse(Int, match(rx, e)[:p])
      (flag, valstring) = Oscar.od_from_p_subgroup(chi, p)
      @test flag
      @test valstring == entry[:valuestring]
    end
  end
end

@testset "order_omega_mod_N" begin
  @test Oscar.order_omega_mod_N(4, 2, 168) == (false, false)
  @test Oscar.order_omega_mod_N(6, 2, 168) == (true, false)
  @test Oscar.order_omega_mod_N(4, 2, 10) == (false, true)
  @test Oscar.order_omega_mod_N(4, 2, 12) == (true, true)
  for (d, q) in [(4, 2), (6, 3), (12, 5), (14, 7)]
    @test Oscar.order_omega_mod_N(d, q, 60)[2]
    @test Oscar.order_omega_mod_N(d, q, ZZ(60))[2]
  end
end

@testset "reduce_mod_squares" begin
  F, z = cyclotomic_field(7)
  for val in [F(0), F(1), -F(5), 1 + 2*z + z^5]
    @test Oscar.reduce_mod_squares(val) == val
    @test Oscar.reduce_mod_squares(12 * val) == 3 * val
    @test Oscar.reduce_mod_squares(27 * val) == 3 * val
    @test Oscar.reduce_mod_squares(val // 12) == 3 * val
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
  @test Oscar._orthogonal_discriminant_indicator0(t[2]) == "-7"
  @test Oscar._orthogonal_discriminant_indicator0(mod(t, 3)[2]) == "O-"
  @test_throws ArgumentError Oscar._orthogonal_discriminant_indicator0(t[4])
end

@testset "orbits under automorphisms" begin
  # _row_orbits_from_column_orbits
  M = identity_matrix(QQ, 5)
  M = M = [[M[i,j] for j in 1:5] for i in 1:5]
  orbs = [[1], [2, 3], [4, 5]]
  @test Oscar._row_orbits_from_column_orbits(M, orbs) == orbs
  M = M[[2, 1, 4, 5, 3]]
  @test Oscar._row_orbits_from_column_orbits(M, orbs) ==        [[1, 5], [2], [3, 4]]

  # _common_orbits
  @test Oscar._common_orbits(
        [[1, 2], [3, 4]], [[1], [2, 3], [4]]) == [[1, 2, 3, 4]]
  @test Oscar._common_orbits(
        [[1], [2, 3], [4]], [[1, 2], [3, 4]]) == [[1, 2, 3, 4]]
  @test Oscar._common_orbits(
        [[1], [2], [3, 4], [5, 6]], [[1, 3], [2, 6], [4], [5]]) ==
        [[1, 3, 4], [2, 5, 6]]
  @test Oscar._common_orbits(
        [[1, 3], [2, 6], [4], [5]], [[1], [2], [3, 4], [5, 6]]) ==
        [[1, 3, 4], [2, 5, 6]]

  # _inverse_fusion
  @test Oscar._inverse_fusion(Int[], 1) == [[]]
  @test Oscar._inverse_fusion([1, 2, 3, 3], 4) ==
        [[1], [2], [3, 4], []]
  @test Oscar._inverse_fusion([1, 3, 3], 4) ==
        [[1], [], [2, 3], []]

  # orbits_group_automorphisms for ordinary tables
  @test Oscar.orbits_group_automorphisms(
        character_table("A5")) == [[1], [2, 3], [4], [5]]
  @test Oscar.orbits_group_automorphisms(
        character_table("A6")) == [[1], [2, 3], [4, 5], [6], [7]]
  @test Oscar.orbits_group_automorphisms(
        character_table("L3(8)")) == [[1], [2], [3, 4, 5, 6, 7, 8],
        [9, 10, 11, 12, 13, 14], [15, 16, 17, 18, 19, 20],
        [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, 32],
        [33], [34, 35, 36], [37, 38, 39, 40, 41, 42],
        [43, 44, 45, 46, 47, 48], [49, 50, 51, 52, 53, 54],
        [55, 56, 57, 58, 59, 60], [61], [62, 63, 64, 65, 66, 67],
        [68, 69], [70, 71, 72]]
  @test Oscar.orbits_group_automorphisms(
        character_table("M11")) == [[x] for x in 1:10]
  @test Oscar.orbits_group_automorphisms(
        character_table("U3(5)")) == [[1], [2], [3], [4, 5, 6], [7], [8],
        [9], [10], [11, 12], [13, 14]]

  # orbits_group_automorphisms for Brauer tables
  @test Oscar.orbits_group_automorphisms(
        character_table("A5", 2)) == [[1], [2, 3], [4]]
  @test Oscar.orbits_group_automorphisms(
        character_table("A6", 2)) == [[1], [2, 3], [4, 5]]
  @test Oscar.orbits_group_automorphisms(
        character_table("L3(8)", 2)) == [[1], [2, 3, 4, 5, 6, 7],
        [8, 9, 10], [11, 12, 13, 14, 15, 16], [17, 18, 19, 20, 21, 22],
        [23, 24, 25, 26, 27, 28], [29, 30, 31, 32, 33, 34],
        [35, 36], [37, 38, 39, 40, 41, 42], [43, 44, 45],
        [46, 47, 48, 49, 50, 51], [52, 53, 54, 55, 56, 57],
        [58, 59, 60, 61, 62, 63], [64]]
  @test Oscar.orbits_group_automorphisms(
        character_table("M11", 2)) == [[x] for x in 1:5]
  @test Oscar.orbits_group_automorphisms(
        character_table("U3(5)", 2)) == [[1], [2], [3, 4, 5], [6],
        [7, 8]]
end

@testset "orthogonal discriminant indicator 0 via a quadratic subfield" begin
  T = character_table(cyclic_group(15))
  chi = first(x for x in T if is_faithful(x))
  # both are correct as sqrt(5) lies in the character field of chi + conj(chi);
  # which one is found depends on the Julia version
  @test Oscar.OrthogonalDiscriminants._orthogonal_discriminant_indicator0(chi) in ("-3", "-15")
end
