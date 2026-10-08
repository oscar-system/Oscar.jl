# For each library, the groups, identifiers and number of groups of a
# selection must agree with each other and with the classic `all_*` function.
function _test_group_library_selection(L, classic, filters)
  S = find(L, filters...)
  ids = collect(keys(S))
  groups = collect(S)

  @test groups isa Vector{eltype(S)}
  @test length(S) == length(ids) == length(groups)
  @test issorted(ids)
  @test [identify(L, G) for G in groups] == ids
  @test [identify(L, G) for G in S] == ids
  @test sort!([identify(L, G) for G in classic(filters...)]) == ids
  @test all(id -> identify(L, L[id]) == id, ids)

  @test isempty(S) == isempty(ids)
  if isempty(ids)
    @test_throws ArgumentError first(S)
  else
    @test identify(L, first(S)) == first(ids)
  end
end

@testset "Group library interface" begin
  @testset "small groups" begin
    L = small_groups_library()
    @testset for filters in [
      (8,),
      (order => 16, !is_abelian),
      (16, exponent => [2, 4]),
      (16, exponent => 5),
      (order => [6, 8, 12], is_abelian),
      (60, is_simple),
      (1:30, number_of_conjugacy_classes => 5),
    ]
      _test_group_library_selection(L, all_small_groups, filters)
    end

    @test describe(L[8, 3]) == "D8"
    @test describe(L[(8, 3)]) == "D8"
    @test identify(L, symmetric_group(4)) == (24, 12)
    @test length(find(L, [8, 16], order => 16:32)) == 14
    @test isempty(find(L, -3:0))

    # counted from the library's tables, there are too many groups to list
    @test length(find(L, order => 512, !is_abelian)) == 10494183
    @test length(find(L, order => 512, is_abelian)) == 30
    @test identify(L, first(find(L, order => 512, !is_abelian))) == (512, 2)
    @test [identify(L, G) for G in Iterators.take(find(L, 512), 3)] ==
      [(512, 1), (512, 2), (512, 3)]

    # the number of groups of order 1024 is known, the groups are not
    @test has_number_of_groups(L, 1024)
    @test !has_groups(L, 1024)
    S = find(L, 1024)
    @test length(S) == 49487367289
    @test_throws ArgumentError first(S)
    @test_throws ArgumentError collect(S)
    @test_throws ArgumentError length(find(L, 1024, is_abelian))

    P = small_groups_library(PermGroup)
    @test P[8, 3] isa PermGroup
    @test collect(find(P, 6)) isa Vector{PermGroup}
    @test first(find(P, 6, !is_abelian)) isa PermGroup

    @test_throws ArgumentError find(L)
    @test_throws ArgumentError find(L, is_abelian)
    @test_throws ArgumentError find(L, 8, gens)
    @test_throws ArgumentError find(L, 8, 5)
    @test_throws ArgumentError L[1, 2]
  end

  @testset "transitive groups" begin
    L = transitive_groups_library()
    @testset for filters in [
      (4,),
      (1,),
      (1, !is_abelian),
      (degree => 1:5, is_abelian),
      (degree => 6, order => 1:12),
      (6, is_primitive),
      (number_of_moved_points => 3:2:9, is_cyclic),
      (8, transitivity => 2),
    ]
      _test_group_library_selection(L, all_transitive_groups, filters)
    end

    @test L[5, 4] == alternating_group(5)
    @test identify(L, symmetric_group(1)) == (1, 1)
    @test length(find(L, 30)) == 5712
    @test collect(keys(find(L, degree => 1:2, is_abelian))) == [(1, 1), (2, 1)]

    @test !has_groups(L, 64)
    S = find(L, 64, is_abelian)
    @test_throws ArgumentError first(S)
    @test_throws ArgumentError length(S)
    @test_throws ArgumentError find(L, is_abelian)
  end

  @testset "primitive groups" begin
    L = primitive_groups_library()
    @testset for filters in [
      (5,),
      (1,),
      (3:2:9, is_abelian),
      (degree => 10, !is_solvable, order => 1:1000),
      (8, transitivity => 3),
      (degree => 2:12, is_cyclic),
    ]
      _test_group_library_selection(L, all_primitive_groups, filters)
    end

    @test order(L[10, 1]) == 60
    @test identify(L, symmetric_group(4)) == (4, 2)
    @test length(find(L, 50)) == 9
  end

  @testset "perfect groups" begin
    L = perfect_groups_library()
    @testset for filters in [
      (order => 1:200,),
      (1:200, !is_simple),
      (1:400, number_of_conjugacy_classes => 5:8),
      (order => 1:5:200, order => 25:60),
      (7200,),
      (17,),
    ]
      _test_group_library_selection(L, all_perfect_groups, filters)
    end

    @test L[120, 1] isa PermGroup
    F = perfect_groups_library(FPGroup)
    @test F[120, 1] isa FPGroup
    @test collect(find(F, 60:120)) isa Vector{FPGroup}
    @test_throws MethodError perfect_groups_library(PcGroup)
    @test_throws ArgumentError find(L, 60, 5)
  end

  @testset "groups with given class number" begin
    L = groups_with_class_number_library()
    @testset for filters in [
      (4,),
      (6, order => 1:20),
      (8, is_abelian),
      (1:3, !is_abelian),
      (number_of_conjugacy_classes => 5, is_simple),
    ]
      _test_group_library_selection(L, all_groups_with_class_number, filters)
    end

    @test number_of_conjugacy_classes(L[5, 6]) == 5
    @test length(find(L, 14)) == number_of_groups_with_class_number(14)
    P = groups_with_class_number_library(PermGroup)
    @test collect(find(P, 3)) isa Vector{PermGroup}
    @test first(find(P, 3, !is_abelian)) isa PermGroup
  end

  # only the GAP package SOTGrps provides the groups of order 2662
  @testset "small groups from SOTGrps" begin
    L = small_groups_library()
    @test has_groups(L, 2662)
    S = find(L, 2662)
    @test length(S) == length(collect(keys(S))) == number_of_small_groups(2662)
    @test order(first(S)) == 2662
    @test identify(L, L[2662, 3]) == (2662, 3)
    @test length(find(L, 2662, is_abelian)) == 3
  end
end
