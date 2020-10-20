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

end
