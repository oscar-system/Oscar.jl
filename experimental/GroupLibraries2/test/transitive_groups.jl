@testset "Transitive groups" begin
   @test has(TransitiveGroupsLibrary, 10)
   @test !has(TransitiveGroupsLibrary, 60)

   @test has_count(TransitiveGroupsLibrary, 10)
   @test !has_count(TransitiveGroupsLibrary, 60)

   @test has_identification(TransitiveGroupsLibrary, 10)
   @test !has_identification(TransitiveGroupsLibrary, 60)

   @test count(TransitiveGroupsLibrary, 4) == 5
   @test count(TransitiveGroupsLibrary, 10) == 45
   @test_throws ArgumentError count(TransitiveGroupsLibrary, 60)

   for i in 1:10
     @test count(TransitiveGroupsLibrary, i) == length(get_all(TransitiveGroupsLibrary, degree => i))
   end

   @test length(get_all(TransitiveGroupsLibrary, degree => -1)) == 0
   @test get_one(TransitiveGroupsLibrary, degree => -1) == nothing

   G = symmetric_group(4)
   H1 = alternating_group(4)
   H2 = sub(G, [G([2, 3, 4, 1])])[1]                   # cyclic
   H3 = sub(G, [G([2, 3, 4, 1]), G([2, 1, 4, 3])])[1]  # dihedral
   H4 = sub(G, [G([3, 4, 1, 2]), G([2, 1, 4, 3])])[1]  # Klein subgroup
   L = [G, H1, H2, H3, H4]
   grps = get_all(TransitiveGroupsLibrary, degree => 4)
   @test get_one(TransitiveGroupsLibrary, degree => 4) in grps
   for K in grps
     @test count(l -> is_isomorphic(K, l), L) == 1
   end
   ids = [identification(TransitiveGroupsLibrary, l) for l in L]
   @test sort(ids) == [(4, i) for i in 1:5]
   @test L == [get(TransitiveGroupsLibrary, id...) for id in ids]
   @test Set(L) == Set(get_all(TransitiveGroupsLibrary, degree => 4))
   @test [H2] == get_all(TransitiveGroupsLibrary, degree => 4, is_cyclic)
   @test [H4] == get_all(TransitiveGroupsLibrary, degree => 4, !is_cyclic, is_abelian)

   grps = get_all(TransitiveGroupsLibrary, degree => 6)
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
     @test length(get_all(TransitiveGroupsLibrary, degree => 6, prop => true)) == count(prop, grps)
     @test length(get_all(TransitiveGroupsLibrary, degree => 6, prop => false)) == count(!prop, grps)

     @test length(get_all(TransitiveGroupsLibrary, degree => 6, prop)) == count(prop, grps)
     @test length(get_all(TransitiveGroupsLibrary, degree => 6, !prop)) == count(!prop, grps)
   end

   @test length(get_all(TransitiveGroupsLibrary, degree => 6, order => 1:12)) == count(g -> order(g) in 1:12, grps)

   @test length(get_all(TransitiveGroupsLibrary, degree => 1:6)) == sum([length(get_all(TransitiveGroupsLibrary, degree => i)) for i in 1:6])

   @test issetequal(get_all(TransitiveGroupsLibrary, 3:2:9), get_all(TransitiveGroupsLibrary, degree => 3:2:9))
   @test issetequal(get_all(TransitiveGroupsLibrary, collect(3:2:9)), get_all(TransitiveGroupsLibrary, 3:2:9))
   @test issetequal(reduce(vcat, (get_all(TransitiveGroupsLibrary, i) for i in 3:2:9)), get_all(TransitiveGroupsLibrary, 3:2:9))

   @test issetequal(get_all(TransitiveGroupsLibrary, 3:2:9, is_abelian), get_all(TransitiveGroupsLibrary, degree => 3:2:9, is_abelian))
   @test issetequal(get_all(TransitiveGroupsLibrary, 9, is_abelian), get_all(TransitiveGroupsLibrary, degree => 9, is_abelian))

   @test_throws ArgumentError get(TransitiveGroupsLibrary, 1, 2)
end
