@testset "Subgroups" begin
   G = symmetric_group(7)
   x = cperm(G,[1,2,3,4,5,6,7])
   y = cperm(G,[1,2,3])
   z = cperm(G,[1,2])

   H,f=sub(G,[x,y])
   K,g=sub(x,y)
   @test H isa PermGroup
   @test H==alternating_group(7)
   @test domain(f)==H
   @test codomain(f)==G
   @test [f(x) for x in gens(H)]==gens(H)
   @test (H,f)==(K,g)
   @test isnormal(G,H)
   H,f=sub(G,[x,z])
   @test H==G
   @test f==id_hom(G)


   G = symmetric_group(4)
   L = subgroups(G)
   @test length(L)==30
   @test L[1] isa Tuple{PermGroup,Oscar.GAPGroupHomomorphism{PermGroup,PermGroup}}
   L1 = [x for x in L if isnormal(G,x[1])]
   K = [H[1] for H in normal_subgroups(G)]
   @test length(K)==4
   for H in L1
      @test H[1] in K
   end
end

@testset "Cosets" begin
  
   G = dihedral_group(8)
   H, mH = center(G)
  
   @test index(G, H) == 4
  
   C = right_coset(H, G[1])
   @test order(C) == length(elements(C))
  
   @test length(right_cosets(G, H)) == index(G, H)
  
   @test length(right_transversal(G, H)) == index(G, H)

   @testset "set comparation for cosets in PermGroup" begin
      G=symmetric_group(5)
      x = G(cperm([1,2,3]))
      y = G(cperm([1,4,5]))
      z = G(cperm([1,2],[3,4]))
      H = sub(G,[y])[1]
      K = sub(G,[z])[1]

      @test_throws ArgumentError rc = right_coset(H,K(z))
      @test_throws ArgumentError rc = left_coset(H,K(z))
      @test_throws ArgumentError dc = double_coset(H,K(z),K)
      @test_throws ArgumentError dc = double_coset(K,K(z),H)
      rc = right_coset(H,x)
      lc = left_coset(H,x)
      dc = double_coset(H,x,K)
      @test acting_domain(rc) == H
      @test representative(rc) == x
      @test representative(lc) == x
      @test Set(elements(rc)) == Set([z for z in rc])          # test iterator
      @test Set(elements(dc)) == Set([z for z in dc])
      @test Set(elements(rc)) == Set([h*x for h in H])
      @test Set(elements(lc)) == Set([x*h for h in H])
      @test Set(elements(dc)) == Set([h*x*k for h in H for k in K])
      @test order(rc) == 3
      @test order(dc) == 6
      r1=rand(rc)
      r2=rand(dc)
      @test r1 in rc
      @test r2 in dc
      @test r1 in dc
      @test issubset(rc,dc)
      @test issubset(left_coset(K,x),dc)
      @test !isbicoset(rc)
   end

end

