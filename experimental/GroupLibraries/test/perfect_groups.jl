@testset "Perfect groups" begin
   G = alternating_group(5)
   @test is_perfect(G)
   @test !is_perfect(symmetric_group(5))

   @test PerfectGroups.get(120, 1) isa PermGroup
   @test PerfectGroups.get(PermGroup, 120, 1) isa PermGroup
   @test PerfectGroups.get(FPGroup, 120, 1) isa FPGroup
   @test_throws ArgumentError PerfectGroups.get(MatGroup, 120, 1)

   @test_throws ArgumentError PerfectGroups.get(17, 0)
   @test_throws ArgumentError PerfectGroups.get(17, 1)
   @test_throws ArgumentError PerfectGroups.get(60, 0)
   @test_throws ArgumentError PerfectGroups.get(60, 2)

   @test is_isomorphic(PerfectGroups.get(60, 1), G)
   @test [PerfectGroups.count(i) for i in 2:59] == [0 for i in 1:58]
   x = PerfectGroups.identification(alternating_group(5))
   @test is_isomorphic(PerfectGroups.get(x[1], x[2]), alternating_group(5))
   @test_throws ArgumentError PerfectGroups.identification(symmetric_group(5))

   @test sum(PerfectGroups.count, 1:59) == 1
   @test PerfectGroups.count(ZZRingElem(60)^3) == 1
   @test_throws ArgumentError PerfectGroups.count(0) # invalid argument
   @test_throws ArgumentError PerfectGroups.count(ZZRingElem(60)^10)  # result not known

   Gs = PerfectGroups.get_all(order => 1:200)
   @test length(Gs) == sum(PerfectGroups.count, 1:200)
   @test Gs == PerfectGroups.get_all(1:200)
   @test length(PerfectGroups.get_all(7200)) == PerfectGroups.count(7200)

   # all_perfect_groups with additional attributes
   @test filter(G -> number_of_conjugacy_classes(G) in 5:8, Gs) == PerfectGroups.get_all(1:200, number_of_conjugacy_classes => 5:8)
   @test filter(is_simple, Gs) == PerfectGroups.get_all(1:200, is_simple)
   @test filter(is_simple, Gs) == PerfectGroups.get_all(1:200, is_simple => true)
   @test filter(!is_simple, Gs) == PerfectGroups.get_all(1:200, !is_simple)
   @test filter(!is_simple, Gs) == PerfectGroups.get_all(1:200, is_simple => false)

   # all_perfect_groups with multiple order specifications
   @test PerfectGroups.get_all(order => 1:5:200, order => 25:50) == PerfectGroups.get_all(order => intersect(1:5:200, 25:50))

   # lazy artifact loading
   @test PerfectGroups.get(1376256, 1) isa PermGroup
end
