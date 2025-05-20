@testset "Riese Schmid groups" begin
   # type Q8
   G = quaternion_group(8)
   @test small_group_identification(G) == (8, 4)
   grps = filter(x -> Oscar._Riese_Schmid_type(x)[1] != "",
                 all_small_groups(1:30))
   @test length(grps) == 1
   @test Oscar._Riese_Schmid_type(grps[1]) == ("Q8", 1)
   @test small_group_identification(grps[1]) == (8, 4)

   # type (Q8, r) or (QD, r)
   @testset for p in [3, 5, 7, 11, 13]
     o = modord(2, p)
     twopart = 1
     while is_even(o)
       o = div(o, 2)
       twopart = twopart * 2
     end
     n = p * 8 * twopart
     ids = map(Oscar._Riese_Schmid_type, all_small_groups(n))
     filt = filter(x -> x[1] != "", ids)
     if p == 7
       # the order of 2 mod 7 is odd
       @test length(filt) == 0
     else
       @test length(filter(x -> x[1] == "Q8", filt)) == 2
       @test length(filter(x -> x[1] == "QD", filt)) > 2
     end
   end
end

@testset "Yamada's examples" begin
   G, chi = Oscar.yamada_example(7, 3)
   @test Oscar.local_schur_indices(chi) == [7 => 3]
   @test schur_index(chi) == 3
   G, chi = Oscar.yamada_example(7, 6)
   @test Oscar.local_schur_indices(chi) == [7 => 6]
   @test schur_index(chi) == 6
   G, chi = Oscar.yamada_example(5, 2)
   @test Oscar.local_schur_indices(chi) == [0 => 2, 5 => 2]
   @test schur_index(chi) == 2

   @test_throws ArgumentError Oscar.yamada_example(2, 2)
   @test_throws ArgumentError Oscar.yamada_example(3, 4)
end

@testset "small examples from Feit's list" begin
   expls = [(190080,  32, [0 => 2, 2 => 2]),     # 2.M12
            (604800, 336, [2 => 2, 3 => 2])]     # J2
   for (grouporder, deg, res) in expls
     G = perfect_group(grouporder, 1)
     t = character_table(G)
     irr = collect(t)
     pos = findfirst(x -> degree(x) == deg, irr)
     chi = irr[pos]
     @test Oscar.local_schur_indices(chi) == res
   end
end
