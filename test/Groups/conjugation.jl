@testset "Conjugacy classes in symmetric groups" begin
  n = 5
  G = symmetric_group(n)
  
  cc = conjugacy_class(G, G[1])
  cc1 = conjugacy_class(G, G[1]^G[2])
  @test order(cc) == order(cc1)
  @test cc == cc1
  
  ccid = conjugacy_class(G,one(G))
  @test order(ccid)==1
  @test elements(ccid) == [one(G)]
  
  x = perm(G,vcat(2:n,[1]))
  cc = conjugacy_class(G,x)
  @test order(cc) == factorial(n-1)
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

  @test order(cc) == 3
  @test Set(elements(cc)) == Set([x^y for y in G])
  y = rand(cc)
  @test y in elements(cc)
  @test order(y) == 2

  C = conjugacy_classes(G)
  @test length(C) == 5
  @test cc in C
  @test sum([order(c) for c in C]) == order(G)
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
end
