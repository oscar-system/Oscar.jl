@testset "Transitive groups" begin
   @test TransitiveGroups.has(10)
   @test !TransitiveGroups.has(60)

   @test TransitiveGroups.has_count(10)
   @test !TransitiveGroups.has_count(60)

   @test TransitiveGroups.has_identification(10)
   @test !TransitiveGroups.has_identification(60)

   @test TransitiveGroups.count(4) == 5
   @test TransitiveGroups.count(10) == 45
   @test_throws ArgumentError TransitiveGroups.count(60)

   for i in 1:10
     @test TransitiveGroups.count(i) == length(TransitiveGroups.get_all(degree => i))
   end

   @test length(TransitiveGroups.get_all(degree => -1)) == 0
   @test TransitiveGroups.get_one(degree => -1) == nothing

   G = symmetric_group(4)
   H1 = alternating_group(4)
   H2 = sub(G, [G([2, 3, 4, 1])])[1]                   # cyclic
   H3 = sub(G, [G([2, 3, 4, 1]), G([2, 1, 4, 3])])[1]  # dihedral
   H4 = sub(G, [G([3, 4, 1, 2]), G([2, 1, 4, 3])])[1]  # Klein subgroup
   L = [G, H1, H2, H3, H4]
   grps = TransitiveGroups.get_all(degree => 4)
   @test TransitiveGroups.get_one(degree => 4) in grps
   for K in grps
     @test count(l -> is_isomorphic(K, l), L) == 1
   end
   ids = [TransitiveGroups.identification(l) for l in L]
   @test sort(ids) == [(4, i) for i in 1:5]
   @test L == [TransitiveGroups.get(id...) for id in ids]
   @test Set(L) == Set(TransitiveGroups.get_all(degree => 4))
   @test [H2] == TransitiveGroups.get_all(degree => 4, is_cyclic)
   @test [H4] == TransitiveGroups.get_all(degree => 4, !is_cyclic, is_abelian)

   grps = TransitiveGroups.get_all(degree => 6)
   @test length(grps) == 16
   props = [
      is_abelian,
      is_almost_simple,
      is_cyclic,
      is_nilpotent,
      is_perfect,
      is_quasisimple,
      is_simple,
      is_sporadic_simple,
      is_solvable,
      is_supersolvable,
      is_transitive,
      is_primitive,
   ]
   @testset "get_all filtering for $(prop)" for prop in props
     @test length(TransitiveGroups.get_all(degree => 6, prop => true)) == count(prop, grps)
     @test length(TransitiveGroups.get_all(degree => 6, prop => false)) == count(!prop, grps)

     @test length(TransitiveGroups.get_all(degree => 6, prop)) == count(prop, grps)
     @test length(TransitiveGroups.get_all(degree => 6, !prop)) == count(!prop, grps)
   end

   @test length(TransitiveGroups.get_all(degree => 6, order => 1:12)) == count(g -> order(g) in 1:12, grps)

   @test length(TransitiveGroups.get_all(degree => 1:6)) == sum([length(TransitiveGroups.get_all(degree => i)) for i in 1:6])

   @test issetequal(TransitiveGroups.get_all(3:2:9), TransitiveGroups.get_all(degree => 3:2:9))
   @test issetequal(TransitiveGroups.get_all(collect(3:2:9)), TransitiveGroups.get_all(3:2:9))
   @test issetequal(reduce(vcat, (TransitiveGroups.get_all(i) for i in 3:2:9)), TransitiveGroups.get_all(3:2:9))

   @test issetequal(TransitiveGroups.get_all(3:2:9, is_abelian), TransitiveGroups.get_all(degree => 3:2:9, is_abelian))
   @test issetequal(TransitiveGroups.get_all(9, is_abelian), TransitiveGroups.get_all(degree => 9, is_abelian))

   @test_throws ArgumentError TransitiveGroups.get(1, 2)
end
