@testset "Transitive groups" begin
   @test number_transitive_groups(4)==5
   @test number_transitive_groups(10)==45

   for i in 1:10
       @test number_transitive_groups(i) == length(all_transitive_groups(degree => i))
   end

   G = symmetric_group(4)
   H1 = alternating_group(4)
   H2 = sub(G,[G([2,3,4,1])])[1]                # cyclic
   H3 = sub(G,[G([2,3,4,1]), G([2,1,4,3])])[1]  # dihedral
   H4 = sub(G,[G([3,4,1,2]), G([2,1,4,3])])[1]  # Klein subgroup
   L = [G,H1,H2,H3,H4]
   grps = all_transitive_groups(degree => 4)
   for K in grps
      @test count(l -> is_isomorphic(K,l), L) == 1
   end
   @test sort([transitive_group_identification(l) for l in L]) == [(4,i) for i in 1:5]
   @test Set(L) == Set(all_transitive_groups(degree => 4))
   @test [H2] == all_transitive_groups(degree => 4, is_cyclic)
   @test [H4] == all_transitive_groups(degree => 4, !is_cyclic, is_abelian)

   grps = all_transitive_groups(degree => 6)
   @test length(grps)==16
   props = [
      is_abelian,
      is_almostsimple,
      is_cyclic,
      is_nilpotent,
      is_perfect,
      is_simple,
      is_solvable,
      is_supersolvable,
      is_transitive,
      is_primitive,
   ]
   @testset "all_transitive_groups filtering for $(prop)" for prop in props
      @test length(all_transitive_groups(degree => 6, prop => true)) == count(prop, grps)
      @test length(all_transitive_groups(degree => 6, prop => false)) == count(!prop, grps)

      @test length(all_transitive_groups(degree => 6, prop)) == count(prop, grps)
      @test length(all_transitive_groups(degree => 6, !prop)) == count(!prop, grps)
   end

   @test length(all_transitive_groups(degree => 6, order => 1:12)) == count(g -> order(g) in 1:12, grps)

   @test length(all_transitive_groups(degree => 1:6)) == sum([length(all_transitive_groups(degree => i)) for i in 1:6])

end

@testset "Transitivity" begin

   G = symmetric_group(4)
   H1 = alternating_group(4)
   H2 = sub(G,[G([2,3,4,1])])[1]                # cyclic
   H3 = sub(G,[G([2,3,4,1]), G([2,1,4,3])])[1]  # dihedral
   H4 = sub(G,[G([3,4,1,2]), G([2,1,4,3])])[1]  # Klein subgroup
   L = [G,H1,H2,H3,H4]

   # FIXME: the following two tests are fishy. The first couple calls
   # really should result in errors ?!?
   @test [transitivity(G,1:i) for i in 1:5]==[1,2,3,4,0]
   @test [transitivity(L[5],1:i) for i in 1:5]==[1,2,1,1,0]

   @test is_transitive(G)
   H = sub(G,[G([2,3,1,4])])[1]
   @test !is_transitive(H)
   @test is_transitive(H,1:3)

   @test [is_semiregular(l) for l in L]==[0,0,1,0,1]
   @test [is_regular(l) for l in L]==[0,0,1,0,1]

   H = sub(G,[G([2,1,4,3])])[1]
   @test is_semiregular(H)
   @test !is_regular(H)
   @test is_regular(H,[1,2])

   @test_throws ArgumentError transitive_group(1, 2)
end

@testset "Perfect groups" begin
   G = alternating_group(5)
   @test is_perfect(G)
   @test !is_perfect(symmetric_group(5))

   @test perfect_group(120,1) isa PermGroup
   @test perfect_group(PermGroup,120,1) isa PermGroup
   @test perfect_group(FPGroup,120,1) isa FPGroup
   @test_throws ArgumentError perfect_group(MatrixGroup,120,1)

   @test_throws ArgumentError perfect_group(17, 0)
   @test_throws ArgumentError perfect_group(17, 1)
   @test_throws ArgumentError perfect_group(60, 0)
   @test_throws ArgumentError perfect_group(60, 2)

   @test is_isomorphic(perfect_group(60,1),G)
   @test [number_perfect_groups(i) for i in 2:59]==[0 for i in 1:58]
   x = perfect_group_identification(alternating_group(5))
   @test is_isomorphic(perfect_group(x[1],x[2]),alternating_group(5))
   @test_throws ErrorException perfect_group_identification(symmetric_group(5))

   @test sum(number_perfect_groups, 1:59) == 1
   @test number_perfect_groups(fmpz(60)^3) == 1
   @test_throws ArgumentError number_perfect_groups(0) # invalid argument
   @test_throws ErrorException number_perfect_groups(fmpz(60)^10)  # result not known
end

@testset "Small groups" begin
   L = all_small_groups(8)
   LG = [abelian_group(PcGroup,[2,4]), abelian_group(PcGroup,[2,2,2]), cyclic_group(8), quaternion_group(8), dihedral_group(8)]
   @test length(L)==5
   @testset for G in LG
      arr = [i for i in 1:5 if is_isomorphic(L[i],G)]
      @test length(arr)==1
      @test small_group_identification(G)==(8,arr[1])
   end
   @test length(all_small_groups(16))==14
   @test length(all_small_groups(order => 16))==14
   @test length(all_small_groups(16, is_abelian))==5
   @test length(all_small_groups(order => 16, !is_abelian))==9
   @test number_small_groups(16)==14
   @test number_small_groups(17)==1

   @test_throws ArgumentError small_group(1, 2)
end

@testset "Primitive groups" begin
   @test has_primitive_groups(50)
   @test_throws ArgumentError primitive_group(1, 1)
   @test number_primitive_groups(50) == 9
end

@testset "Atlas groups" begin
   # `atlas_group` for type and group name
   @test order(atlas_group(PermGroup, "A5")) == 60
   @test order(atlas_group(MatrixGroup, "A5")) == 60
   @test_throws ErrorException atlas_group(PermGroup, "B")

   # prescribe permutation degree
   info1 = all_atlas_group_infos("A5", degree => 5)
   info2 = all_atlas_group_infos("A5", degree => [5, 6])
   @test length(info1) == 1
   @test length(info2) > length(info1)

   # prescribed character of the matrix group
   t = character_table("A5", 2)
   info1 = all_atlas_group_infos("A5", character => t[2])
   @test length(info1) == 1
   @test info1[1][:dim] == 2
#T info2 = all_atlas_group_infos("A5", character => [t[2], t[3]])
#T not yet supported in GAP

   # prescribed characteristic of the matrix group
   info1 = all_atlas_group_infos("A5", characteristic => 2)
   info2 = all_atlas_group_infos("A5", characteristic => [2, 3])
   @test length(info1) > 0
   @test length(info2) > length(info1)

   # prescribed dimension of the matrix group
   info1 = all_atlas_group_infos("A5", dim => 2)
   info2 = all_atlas_group_infos("A5", dim => [2, 3])
   @test length(info1) > 0
   @test length(info2) > length(info1)

   # prescribed `is_transitive` of the permutation action
   info = all_atlas_group_infos("A5", is_transitive)
   @test length(info) > 0
   @test info == all_atlas_group_infos("A5", is_transitive => true)
   info = all_atlas_group_infos("M11xA6.2^2", !is_transitive)
   @test length(info) > 0
   @test info == all_atlas_group_infos("M11xA6.2^2", is_transitive => false)

   # prescribed `is_primitive` of the permutation action
   info = all_atlas_group_infos("A5", is_primitive)
   @test length(info) > 0
   info = all_atlas_group_infos("2.A5", !is_primitive)
   @test length(info) > 0

   # prescribed rank of the permutation action
   info1 = all_atlas_group_infos("A5", rank_action => 2)
   info2 = all_atlas_group_infos("A5", rank_action => [2, 3])
   @test length(info1) > 1
   @test length(info2) > length(info1)
   @test info1[1][:degree] == 5

   # prescribed transitivity of the permutation action
   info1 = all_atlas_group_infos("A5", transitivity => 2)
   info2 = all_atlas_group_infos("A5", transitivity => 2:3)
   @test length(info1) == 1
   @test length(info2) > length(info1)
   @test info2[1][:degree] == 5

   # prescribed `base_ring` of the matrix group
   # (for different types of fields)
   rings = [GF(2), GF(2, 2), GF(fmpz(2), 2)]
   for R in rings
     info = all_atlas_group_infos("A5", base_ring => R)
     @test length(info) > 0
     @test all(x -> domain(x[:base_ring_iso]) == R, info)

     # create a group from the info
     G = atlas_group(info[1])
     @test order(G) == 60
     @test base_ring(G) == R
   end
end
