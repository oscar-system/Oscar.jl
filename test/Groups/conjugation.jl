@testset "Conjugacy classes in symmetric groups" begin
  n = 5
  G = symmetric_group(n)
  
  cc = conjugacy_class(G, G[1])
  cc1 = conjugacy_class(G, G[1]^G[2])
  @test length(cc) == length(cc1)
  @test cc == cc1
  
  ccid = conjugacy_class(G,one(G))
  @test length(ccid)==1
  @test elements(ccid) == [one(G)]
  
  x = perm(G,vcat(2:n,[1]))
  cc = conjugacy_class(G,x)
  @test length(cc) == factorial(n-1)
  @test x == representative(cc)
  y = rand(cc)
  @test order(y) == order(x)
  @test isconjugate(G,x,x^y)[1]
  z = isconjugate(G,x,x^y)[2]
  @test x^z == x^y
  z = cperm(G,[1,2])
  @test !isconjugate(G,x,z)[1]


# something in smaller dimension
  @test number_conjugacy_classes(symmetric_group(4)) == 5
  G = symmetric_group(4)
  x = perm(G,[3,4,1,2])
  cc = conjugacy_class(G,x)

  @test length(cc) == 3
  @test Set(elements(cc)) == Set([x^y for y in G])
  y = rand(cc)
  @test y in elements(cc)
  @test order(y) == 2

  C = conjugacy_classes(G)
  @test length(C) == 5
  @test cc in C
  @test sum([length(c) for c in C]) == order(G)
  @test sum([x in elements(c) for c in C]) == 1          # x belongs to a unique conjugacy class
  @test sum([y in elements(c) for c in C]) == 1          # x belongs to a unique conjugacy class
  z = rand(G)
  @test sum([z in elements(c) for c in C]) == 1          # x belongs to a unique conjugacy class
  @testset for i in 1:5
     c = C[i]
     x = rand(c)
     y = rand(c)
     @test isconjugate(G,x,y)[1]
     z = isconjugate(G,x,y)[2]
     @test x^z == y
     y = rand(C[(i%5)+1])
     @test !isconjugate(G,x,y)[1]
  end

  CC = conjugacy_classes_subgroups(G)
  @test length(CC)==11
  @testset for i in 1:length(CC)
     @test CC[i] == conjugacy_class(G,representative(CC[i]))
     @test length(CC[i]) == index(G, normalizer(G,representative(CC[i]))[1])
  end
  H=rand(subgroups(G))
  @test sum([length(c) for c in CC]) == length(subgroups(G))
  @test sum([H in elements(c) for c in CC]) == 1         # H belongs to a unique conjugacy class
  @testset for i in 1:length(CC)
     c = CC[i]
     x = rand(c)
     y = rand(c)
     @test isconjugate(G,x,y)[1]
     z = isconjugate(G,x,y)[2]
     @test x^z == y
     y = rand(CC[(i % length(CC))+1])
     @test !isconjugate(G,x,y)[1]
  end

  CC = conjugacy_classes_maximal_subgroups(G)
  @test length(CC)==3
  @test Set([order(representative(l)) for l in CC])==Set([6,8,12])

  x = G(cperm([1,2,3,4]))
  H = sub(G,[x])[1]
  @test normalizer(G,H)==normalizer(G,x)

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
   @testset for y in elements(cc)
       @test sum([y==x^t for t in T])==1
   end
   
   cs = conjugacy_class(G,Cx)
   Nx = normalizer(G,Cx)[1]
   @test index(G,Nx)==length(cs)
   T=right_transversal(G,Nx)
   # Set([Cx^t for t in T]) == Set(elements(cs)) does not work
   @testset for H in elements(cs)
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
   R,t = PolynomialRing(base_ring(G),"t")

   x = generalized_jordan_block(t-1,8)
   y = generalized_jordan_block(t-1,8)
   x[7,8]=gen(base_ring(G))
   x=S(x); y=S(y);
   vero, z = isconjugate(G,x,y)
   @test vero
   @test z in G
   @test x^z==y
   vero, z = isconjugate(S,x,y)
   @test !vero
   x.elm[7,8]=gen(base_ring(G))^8
   vero, z = isconjugate(S,x,y)
   @test z in S
   @test x^z==y

   G = GL(8,5)
   S = SL(8,5)
   R,t = PolynomialRing(base_ring(G),"t")
   x = diagonal_join(generalized_jordan_block(t-1,4), generalized_jordan_block(t-1,2), identity_matrix(base_ring(G),2))
   C = centralizer(G,G(x))[1]
   @test order(C) == order(GL(2,5))*4^2*5^16
   Cs = centralizer(S,S(x))[1]
   @test order(Cs) == div(order(GL(2,5))*4^2*5^16,4)
   x = diagonal_join( [generalized_jordan_block(t-1,2) for i in 1:4] )
   C = centralizer(G,G(x))[1]
   @test order(C) == order(GL(4,5))*5^16
   Cs = centralizer(S,S(x))[1]
   @test order(Cs) == div(order(GL(4,5))*5^16,2)
   x = diagonal_join( [generalized_jordan_block(t-1,4) for i in 1:2] )
   C = centralizer(G,G(x))[1]
   @test order(C) == order(GL(2,5))*5^12
   Cs = centralizer(S,S(x))[1]
   @test order(Cs) == order(GL(2,5))*5^12
   

   x = companion_matrix(t^8+t^3+t^2+t+2)
   C = centralizer(G,G(x))[1]
   @test order(C)==5^8-1
end

@testset "Jordan structure" begin
   F = GF(3,1)[1]
   R,t = PolynomialRing(F,"t")

   L_big = [
        [(t-1,3), (t^2+1,1), (t^2+1,2)],
        [(t^9 + t^7 + 2*t^6 + t^5 + 2*t^4 + t^3 + 2*t^2 + 2*t + 1, 1)],
        [(t-2,9)]
       ]
   for L in L_big
      x = diagonal_join([generalized_jordan_block(a...) for a in L])
      G = GL(9,F)
   # TODO: these will work when the polynomial rings are recognized as the same
   #   @test pol_elementary_divisors(x)==L
   #   @test pol_elementary_divisors(G(x))==L
      s, u = multiplicative_jordan_decomposition(G(x))
      @test parent(s)==G
      @test parent(u)==G
      @test iscoprime(order(s),3)
      @test isone(u) || ispower(order(u))[2]==3
      @test issemisimple(s)
      @test isunipotent(u)
      @test s*u==G(x)
      @test s*u==u*s

      z = rand(G).elm
      x = z^-1*x*z
      a,b = generalized_jordan_form(x)
      @test b^-1*a*b==x
      z = rand(G).elm
      @test generalized_jordan_form(z^-1*x*z)[1]==a
   end
   
   x = one(G)
   @test issemisimple(x) && isunipotent(x)
end
