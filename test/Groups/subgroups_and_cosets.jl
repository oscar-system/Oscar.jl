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
   
   @test derived_subgroup(G)[1] == alternating_group(4)
   L = derived_series(G)
   @test L[1][1] == G
   @test L[2][1] == alternating_group(4)
   @test L[3][1] == sub(G, [cperm([1,2],[3,4]), cperm([1,3],[2,4])])[1]
   @test L[4][1] == sub(G, [one(G)])[1]
end

@testset "Centralizers and Normalizers in Sym(n)" begin
   G = symmetric_group(6)
   x = cperm(G,[1,2])
   y = cperm(G,[1,2,3,4])

   Cx = centralizer(G,x)[1]
   Cy = centralizer(G,y)[1]
   @testset for i in 1:ngens(Cx)
      @test Cx[i]*x == x*Cx[i]
   end
   @testset for i in 1:ngens(Cy)
      @test Cy[i]*y == y*Cy[i]
   end
   notx = setdiff(G,Cx)
   noty = setdiff(G,Cy)
   @testset for i in 1:3
      z=rand(notx)
      @test z*x != x*z
      z=rand(noty)
      @test z*y != y*z
   end
   @test x in Cx
   @test one(G) in Cx
   @test isisomorphic(Cx, direct_product(symmetric_group(2),symmetric_group(4))[1])[1]
   @test isisomorphic(Cy, direct_product(sub(G,[y])[1], symmetric_group(2))[1])[1]

   Nx = normalizer(G,Cx)[1]
   Ny = normalizer(G,Cy)[1]
   @test isnormal(Nx,Cx)
   @test isnormal(Ny,Cy)
   notx = setdiff(G,Nx)
   noty = setdiff(G,Ny)
   @testset for i in 1:3
      z=rand(notx)
      @test Cx^z != Cx
      z=rand(noty)
      @test Cy^z != Cy
   end

   CCx = centralizer(G,Cx)[1]
   CCy = centralizer(G,Cy)[1]
   @testset for i in 1:ngens(Cx)
      for j in 1:ngens(CCx)
         @test Cx[i]*CCx[j] == CCx[j]*Cx[i]
      end
   end
   @testset for i in 1:ngens(Cy)
      for j in 1:ngens(CCy)
         @test Cy[i]*CCy[j] == CCy[j]*Cy[i]
      end
   end
   notx = setdiff(G,Nx)
   noty = setdiff(G,Ny)
   @testset for i in 1:3
      z=rand(notx)
      @testset for j in 1:ngens(Cx)
         @test z*Cx[j] != Cx[j]*z
      end
      z=rand(noty)
      @testset for j in 1:ngens(Cy)
         @test z*Cy[j] != Cy[j]*z
      end
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

