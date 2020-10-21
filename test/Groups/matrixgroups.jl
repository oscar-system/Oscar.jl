@testset "Fields assignment" begin
   T,t=PolynomialRing(FiniteField(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   @test G isa MatrixGroup
   @test F==base_ring(G)
   @test 2==degree(G)
   @test !isdefined(G,:X)
   @test !isdefined(G,:gens)

   @test order(G)==5760
   @test isdefined(G,:X)
   @test !isdefined(G,:gens)

   @test ngens(G)==2
   @test isdefined(G,:gens)

   x = G[1]
   @test !isdefined(x,:X)

   x = matrix(F,2,2,[1,0,0,1])
   x = G(x)
   @test !isdefined(x,:X)
   x = G[1].elm
   x = G(x)
   @test !isdefined(x,:X)
   x = matrix(F,2,2,[1,z,z,1])
   x = G(x)
   @test !isdefined(x,:X)
   @test order(x)==8
   @test isdefined(x,:X)
end

@testset "Membership" begin
   T,t=PolynomialRing(FiniteField(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   S = SL(2,F)
   O = GO(1,2,F)
   
   x = matrix(F,2,2,[1,z,0,z])
   @test x in G
   @test !(x in S)
   @test_throws ArgumentError S(x)
   x = G(x)
   @test parent(x)==G
   @test x==G(x)

   x = matrix(F,2,2,[1,z,0,1])
   @test x in S
   x = S(x)
   @test parent(x)==S
   @test x in G
   @test parent(G(x))==G
   @test parent(S(x))==S
   x = (O[1]*O[2]).elm
   @test x in G
   @test x in O
   @test parent(O(x))==O
   x = G(x)
   @test parent(O(x))==O

end
