@testset "Transitivity" begin
   @test number_transitive_groups(10)==45
   
   @test number_transitive_groups(4)==5
   G = symmetric_group(4)
   H1 = alternating_group(4)
   H2 = sub(G,[G([2,3,4,1])])[1]                 #cyclic
   H3 = sub(G,[G([2,3,4,1]), G([2,1,4,3])])[1]         #dihedral
   H4 = sub(G,[G([3,4,1,2]), G([2,1,4,3])])[1]         #Klein subgroup
   L = [G,H1,H2,H3,H4]
   for i in 1:5
      @test sum([isisomorphic(transitive_group(4,i),l) for l in L])==1
   end
   @test Set([transitive_identification(l) for l in L])==Set(1:5)
   @test Set(L)==Set(all_transitive_groups(degree,4))
   @test [H2]==all_transitive_groups(degree,4, iscyclic)
   @test [H4]==all_transitive_groups(degree,4, iscyclic, false, isabelian)
   @test length(all_transitive_groups(degree,6))==16

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

   @test_throws AssertionError transitive_group(1, 2)
end

@testset "Perfect groups" begin
   G = alternating_group(5)
   @test isperfect(G)
   @test !isperfect(symmetric_group(5))

   @test perfect_group(120,1) isa PermGroup
   @test perfect_group(FPGroup,120,1) isa FPGroup
   @test_throws ArgumentError perfect_group(MatrixGroup,120,1)

   @test isisomorphic(perfect_group(60,1),G)
   @test [number_perfect_groups(i) for i in 2:59]==[0 for i in 1:58]
   x = perfect_identification(alternating_group(5))
   @test isisomorphic(perfect_group(x[1],x[2]),alternating_group(5))

   @test_throws AssertionError perfect_group(60, 2)
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
   @test length(all_small_groups(16, isabelian))==5
   @test number_small_groups(16)==14
   @test number_small_groups(17)==1   

   @test_throws AssertionError small_group(1, 2)
end

@testset "Primitive groups" begin
   @test_throws AssertionError primitive_group(1, 1)
end
