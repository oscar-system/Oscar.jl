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
      @test count(l -> isisomorphic(K,l), L) == 1
   end
   @test sort([transitive_group_identification(l) for l in L]) == [(4,i) for i in 1:5]
   @test Set(L) == Set(all_transitive_groups(degree => 4))
   @test [H2] == all_transitive_groups(degree => 4, iscyclic)
   @test [H4] == all_transitive_groups(degree => 4, !iscyclic, isabelian)

   grps = all_transitive_groups(degree => 6)
   @test length(grps)==16
   props = [
      isabelian,
      isalmostsimple,
      iscyclic,
      isnilpotent,
      isperfect,
      issimple,
      issolvable,
      issupersolvable,
      istransitive,
      isprimitive,
   ]
   @testset "all_transitive_groups filtering for $(prop)" for prop in props
      @test length(all_transitive_groups(degree, 6, prop, true)) == count(prop, grps)
      @test length(all_transitive_groups(degree, 6, prop, false)) == count(!prop, grps)

      @test length(all_transitive_groups(degree, 6, prop => true)) == count(prop, grps)
      @test length(all_transitive_groups(degree, 6, prop => false)) == count(!prop, grps)

      @test length(all_transitive_groups(degree, 6, prop)) == count(prop, grps)
      @test length(all_transitive_groups(degree, 6, !prop)) == count(!prop, grps)
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

   @test istransitive(G)
   H = sub(G,[G([2,3,1,4])])[1]
   @test !istransitive(H)
   @test istransitive(H,1:3)

   @test [issemiregular(l) for l in L]==[0,0,1,0,1]
   @test [isregular(l) for l in L]==[0,0,1,0,1]

   H = sub(G,[G([2,1,4,3])])[1]
   @test issemiregular(H)
   @test !isregular(H)
   @test isregular(H,[1,2])

   @test_throws ArgumentError transitive_group(1, 2)
end

@testset "Perfect groups" begin
   G = alternating_group(5)
   @test isperfect(G)
   @test !isperfect(symmetric_group(5))

   @test perfect_group(120,1) isa PermGroup
   @test perfect_group(PermGroup,120,1) isa PermGroup
   @test perfect_group(FPGroup,120,1) isa FPGroup
   @test_throws ArgumentError perfect_group(MatrixGroup,120,1)

   @test_throws ArgumentError perfect_group(17, 0)
   @test_throws ArgumentError perfect_group(17, 1)
   @test_throws ArgumentError perfect_group(60, 0)
   @test_throws ArgumentError perfect_group(60, 2)

   @test isisomorphic(perfect_group(60,1),G)
   @test [number_perfect_groups(i) for i in 2:59]==[0 for i in 1:58]
   x = perfect_group_identification(alternating_group(5))
   @test isisomorphic(perfect_group(x[1],x[2]),alternating_group(5))
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
      arr = [i for i in 1:5 if isisomorphic(L[i],G)]
      @test length(arr)==1
      @test small_group_identification(G)==(8,arr[1])
   end
   @test length(all_small_groups(16))==14
   @test length(all_small_groups(order => 16))==14
   @test length(all_small_groups(16, isabelian))==5
   @test length(all_small_groups(order => 16, !isabelian))==9
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
   @test_throws ErrorException atlas_group(PermGroup, "B")
end
