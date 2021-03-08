@testset "Finitely presented groups" begin
   F = free_group(2)
   x,y = gens(F)
   @test x == F[1]
   @test y == F[2]
   
   n=5
   G,f = quo(F, [x^2,y^n,(x*y)^2] )           # dihedral group D(2n)
   @test isfinite(G)
   @test order(G) == 2*n
   @test !isabelian(G)
   @test isisomorphic(G, dihedral_group(2*n))[1]
   @test !isinjective(f)
   @test issurjective(f)
   @test exponent(G) == 2*n
   if order(G[2])==n
      @test G[2]^G[1] == G[2]^-1
   else
      @test G[1]^G[2] == G[1]^-1
   end
   
   G,f = quo(F, [x^n,y^n,comm(x,y)])          # group C(n) x C(n)
   @test isfinite(G)
   @test order(G) == n^2
   @test isabelian(G)
   @test !isinjective(f)
   @test issurjective(f)
   @test exponent(G) == n
   @test isone(G[1]^n)
   @test relators(G)==[x^n,y^n,comm(x,y)]

   S = symmetric_group(4)
   G,f = quo(S, [cperm(S,[1,3,2])])
   @test order(G) == 2
   @test f(S([1,2,4,3]))==G[1]
   @test f(S([2,1,4,3]))==one(G)
end
