@testset "Permutations" begin
  for n = 10:13
    G=symmetric_group(n)
    x=cperm(1:5,6:8,9:n)
    A=vcat([i for i in 10:n],[9])
    A=vcat([2,3,4,5,1,7,8,6],A)
    y=G(A)

    @test x==y
    @test A==listperm(y)
    @test x==G(listperm(x))
    @test order(x) == lcm(15,n-8)
    @test x(3)==4
    @test x(8)==6
    @test x(n)==9
  end
end
 
@testset "Change of parent" begin
  for n = 10:12
    H = symmetric_group(n+3)
    K = symmetric_group(n-1)
    x = cperm(H, 1:5,6:8,9:n)
    @test parent(x)===H
    @test_throws ArgumentError cperm(K,1:5,6:8,9:n)
    y = H(x)
    @test parent(y) == H
    @test parent(y) === H
    @test parent(x) === H
    @test_throws ArgumentError K(x)

#    z = perm(vcat(2:(n-2),[1]))
#    @test parent(z) == symmetric_group(n-2)
    z = perm(symmetric_group(n),vcat(2:(n-2),[1]))
    @test parent(z) == symmetric_group(n)
    @test_throws MethodError perm(symmetric_group(n-3),z)

    @test cperm(K,Int64[]) == one(K)
  end
end

@testset "Generators" begin
   L=[symmetric_group(4), cyclic_group(5), free_group(3), symplectic_group(4,3)]
   for G in L
      K=gens(G)
      @test length(K) == ngens(G)
      for i in 1:length(K)
         @test K[i] == G[i]
         @test K[i] == gens(G,i)
      end
   end
end
   
