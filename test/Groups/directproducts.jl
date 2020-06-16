@testset "Directproducts" begin
   S = symmetric_group(4)
   C = abelian_group(PcGroup, [2,2])
   G = directproduct(S,C)

   @test G isa DirectProductOfGroups
   @test order(G)==order(S)*order(C)
   @test exponent(G)==lcm(exponent(S),exponent(C))
   @test typeof(rand(G))==Oscar.GAPGroupElem{DirectProductOfGroups{PermGroup,PcGroup}}
   @test factorofdirectproduct(G,1)==S
   @test factorofdirectproduct(G,2)==C

   x = rand(G)
   @test x in G
   @test projection(G,1)(x) in S
   @test projection(G,2)(x) in C
   @test embedding(G,1)(rand(S)) in G
   @test embedding(G,2)(rand(C)) in G
   @test x==G(projection(G,1)(x), projection(G,2)(x))
   @test x==embedding(G,1)(projection(G,1)(x))*embedding(G,2)(projection(G,2)(x))
   S1 = image(embedding(G,1))[1]
   C1 = image(embedding(G,2))[1]
   @test intersect(G,S1)[1]==S1
   @test isisomorphic(S1,S)[1]
   @test isisomorphic(quo(G,S1)[1],C)[1]
   @test isisomorphic(quo(G,C1)[1],S1)[1]

   x = G(cperm([1,2]),C[1])
   H = sub(G,[x])[1]
   @test order(H)==2
   @test index(G,H) isa Integer
   @test_throws ArgumentError write_as_full(H)
   @test intersect(G,H)[1]==H
   x = G(cperm([1,2,3]),C[1])
   H = sub(G,[x])[1]
   @test order(H)==6
   @test H==write_as_full(H)
   @test intersect(G,H)[1]==H

   P1=sylow_subgroup(S,2)[1]
   P2=sylow_subgroup(C,2)[1]
   P=sylow_subgroup(G,2)[1]
   PP=directproduct(P1,P2)
   @test order(P)==order(P1)*order(P2)
   @test isconjugate(G,P,PP)[1]
   x = isconjugate(G,P,PP)[2]
   @test P^x==PP
end
