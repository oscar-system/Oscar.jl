@testset "Conjugacy classes in symmetric groups" begin
  n = 5
  G = symmetric_group(n)
  
  cc = conjugacy_class(G, G[1])
  @test acting_group(cc) === G
  @test representative(cc) == G[1]
  cc1 = conjugacy_class(G, G[1]^G[2])
  @test length(cc) == length(cc1)
  @test length(cc) == length(cc1)
  @test cc == cc1
  
  ccid = conjugacy_class(G,one(G))
  @test acting_group(ccid) === G
  @test representative(ccid) == one(G)
  @test length(ccid)==1
  @test collect(ccid) == [one(G)]
  
  x = perm(G,vcat(2:n,[1]))
  cc = conjugacy_class(G,x)
  @test acting_group(cc) === G
  @test representative(cc) == x
  @test length(cc) == factorial(n-1)
  @test x == representative(cc)
  y = rand(cc)
  @test order(y) == order(x)
  @test isconjugate(G,x,x^y)
  @test representative_action(G,x,x^y)[1]
  z = representative_action(G,x,x^y)[2]
  @test x^z == x^y
  z = cperm(G,[1,2])
  @test !isconjugate(G,x,z)
  @test !representative_action(G,x,z)[1]


  @inferred fmpz number_conjugacy_classes(symmetric_group(4))
  @inferred fmpz number_conjugacy_classes(symmetric_group(40))

# something in smaller dimension
  @test number_conjugacy_classes(symmetric_group(4)) == 5
  G = symmetric_group(4)
  x = perm(G,[3,4,1,2])
  cc = conjugacy_class(G,x)

  @test length(cc) == 3
  @test Set(collect(cc)) == Set(x^y for y in G)
  y = rand(cc)
  @test y in collect(cc)
  @test order(y) == 2

  C = conjugacy_classes(G)
  @test length(C) == 5
  @test all(cc -> acting_group(cc) === G, C)
  @test cc in C
  @test sum(length, C) == order(G)
  @test count(c -> x in c, C) == 1          # x belongs to a unique conjugacy class
  @test count(c -> y in c, C) == 1          # y belongs to a unique conjugacy class
  z = rand(G)
  @test count(c -> z in c, C) == 1          # z belongs to a unique conjugacy class
  @testset for i in 1:5
     c = C[i]
     x = rand(c)
     y = rand(c)
     @test isconjugate(G,x,y)
     @test representative_action(G,x,y)[1]
     z = representative_action(G,x,y)[2]
     @test x^z == y
     y = rand(C[(i%5)+1])
     @test !isconjugate(G,x,y)
     @test !representative_action(G,x,y)[1]
  end

  CC = @inferred conjugacy_classes_subgroups(G)
  @test length(CC)==11
  @test all(cc -> acting_group(cc) === G, CC)
  @testset for C in CC
     @test C == conjugacy_class(G, representative(C))
     @test length(C) == index(G, normalizer(G, representative(C))[1])
     @test degree(representative(C)) == degree(G)
  end
  H=rand(subgroups(G))
  @test sum([length(c) for c in CC]) == length(subgroups(G))
  @test count(c -> H in c, CC) == 1          # H belongs to a unique conjugacy class
  @testset for i in 1:length(CC)
     c = CC[i]
     x = rand(c)
     y = rand(c)
     @test isconjugate(G,x,y)
     @test representative_action(G,x,y)[1]
     z = representative_action(G,x,y)[2]
     @test x^z == y
     y = rand(CC[(i % length(CC))+1])
     @test !isconjugate(G,x,y)
     @test !representative_action(G,x,y)[1]
  end

  CC = @inferred conjugacy_classes_maximal_subgroups(G)
  @test length(CC)==3
  @test Set([order(Int, representative(l)) for l in CC])==Set([6,8,12])

  x = G(cperm([1,2,3,4]))
  H = sub(G,[x])[1]
  @test normalizer(G,H)==normalizer(G,x)

  G = symmetric_group(5)
  CC = @inferred conjugacy_classes_maximal_subgroups(G)
  all(H -> degree(H) == degree(G), map(representative, CC))

  G = symmetric_group(10)
  x = rand(G)
  H = sub(G,[x])[1]
  z = rand(G)
  K = H^z
  @test Set([G(y) for y in K]) == Set([G(y^z) for y in H])
#  @test Set(K) == Set([y^z for y in H])  may not work because the parent of the elements are different

end

function TestConjCentr(G,x)
   Cx = centralizer(G,x)[1]
   cc = conjugacy_class(G,x)
   @test index(G,Cx)==length(cc)
   T=right_transversal(G,Cx)
   @testset for y in cc
       @test count(t -> y==x^t, T) == 1
   end
   
   cs = conjugacy_class(G,Cx)
   Nx = normalizer(G,Cx)[1]
   @test index(G,Nx)==length(cs)
   T=right_transversal(G,Nx)
   # Set([Cx^t for t in T]) == Set(collect(cs)) does not work
   @testset for H in collect(cs)
       @test sum([H==Cx^t for t in T])==1
   end
end

@testset "Conjugation and centralizers" begin
   G = symmetric_group(6)
   x = cperm(G,[1,2,3,4])
   TestConjCentr(G,x)

   G = GL(3,3)
   x = G[2]
   TestConjCentr(G,x)

   G = special_unitary_group(2,3)
   x = G[1]
   TestConjCentr(G,x)
end

@testset "Conjugation and centralizers for GL and SL" begin
   G = GL(8,25)
   S = SL(8,25)
   l = gen(base_ring(G))
   R,t = PolynomialRing(base_ring(G),"t")

   x = generalized_jordan_block(t-1,8)
   y = generalized_jordan_block(t-1,8)
   x[7,8]=l
   x=S(x); y=S(y);
   vero, z = isconjugate(G,x,y)
   @test vero
   @test z in G
   @test x^z==y
   vero, z = isconjugate(S,x,y)
   @test !vero
   x.elm[7,8]=l^8
   vero, z = isconjugate(S,x,y)
   @test z in S
   @test x^z==y

   G = GL(8,5)
   S = SL(8,5)
   R,t = PolynomialRing(base_ring(G),"t")
   x = cat(generalized_jordan_block(t-1,4), generalized_jordan_block(t-1,2), identity_matrix(base_ring(G),2), dims=(1,2))
   C = centralizer(G,G(x))[1]
   @test order(C) == order(GL(2,5))*4^2*5^16
   @testset for y in gens(C)
      @test x*y==y*x
   end
   Cs = centralizer(S,S(x))[1]
   @test order(Cs) == div(order(GL(2,5))*4^2*5^16,4)
   @testset for y in gens(Cs)
      @test x*y==y*x
   end
   x = cat( [generalized_jordan_block(t-1,2) for i in 1:4]..., dims=(1,2) )
   C = centralizer(G,G(x))[1]
   @test order(C) == order(GL(4,5))*5^16
   @testset for y in gens(C)
      @test x*y==y*x
   end
   Cs = centralizer(S,S(x))[1]
   @test order(Cs) == div(order(GL(4,5))*5^16,2)
   x = cat( [generalized_jordan_block(t-1,4) for i in 1:2]..., dims=(1,2) )
   C = centralizer(G,G(x))[1]
   @test order(C) == order(GL(2,5))*5^12
   Cs = centralizer(S,S(x))[1]
   @test order(Cs) == order(GL(2,5))*5^12
   

   x = companion_matrix(t^8+t^3+t^2+t+2)
   C = centralizer(G,G(x))[1]
   @test order(C)==5^8-1
   @testset for y in gens(C)
      @test x*y==y*x
   end

   x = cat(companion_matrix((t^2+3)^2), companion_matrix(t^3+3*t+2), companion_matrix(t-3), dims=(1,2))
   C = centralizer(G,G(x))[1]
   @test order(C)==24*124*4*25
   @testset for y in gens(C)
      @test x*y==y*x
   end
   C = centralizer(S,S(x))[1]
   @test order(C)==24*124*25
   @testset for y in gens(C)
      @test x*y==y*x
   end

   F,t = PolynomialRing(GF(3),"t")
   F,z = FiniteField(t^2+1,"z")
   _,t = PolynomialRing(F,"t")
   G = GL(8,F)
   x = cat(generalized_jordan_block(t^2+t+z,2), generalized_jordan_block(t^2+z+1,2); dims=(1,2))
   C = centralizer(G,G(x))[1]
   @test order(C)==81^2*80^2
   @testset for y in gens(C)
      @test x*y==y*x
   end

end
