@testset "Permutations" begin
  for n = 10:13
    G=symmetric_group(n)
    x=cperm(1:5,6:8,9:n)
    A=vcat([i for i in 10:n],[9])
    A=vcat([2,3,4,5,1,7,8,6],A)
    y=G(A)

    @test x==y
    @test A==Vector(y)
    @test typeof(Vector(y))==Vector{Int64}
    @test typeof(Vector{fmpz}(y))==Vector{fmpz}
    @test x==G(Vector(x))
    @test is_finiteorder(x)
    @test order(x) == lcm(15,n-8)
    for T in [Int, BigInt, fmpz]
      @test order(T, x) == lcm(15,n-8)
      @test order(T, x) isa T
    end

    @test x(3)==4
    @test x(8)==6
    @test x(n)==9

    @test 3^x == 4
    @test n^x == 9

    for T in [Int32, Int, BigInt, fmpz]
      @test x(T(n)) == T(9)
      @test typeof(x(T(n))) == T
      @test T(n)^x == T(9)
      @test typeof(T(n)^x) == T
    end
  end

  G=symmetric_group(6)
  x=gap_perm([2,3,4,5,6,1])
  @test x==gap_perm(Int8[2,3,4,5,6,1])
  @test x==gap_perm(fmpz[2,3,4,5,6,1])
  @test x==perm(G,[2,3,4,5,6,1])
  @test x==perm(G,Int8[2,3,4,5,6,1])
  @test x==perm(G,fmpz[2,3,4,5,6,1])
  @test cperm(G,Int[])==one(G)
  @test x==cperm(G,1:6)
  @test x==cperm(G,[1,2,3,4,5,6])
  @test one(G)==gap_perm(1:6)
  @test_throws ArgumentError G(gap_perm([2,3,4,5,6,7,1]))
  @test_throws ArgumentError G([2,3,1,4,6,5,7])
  @test G(gap_perm([2,3,1,4,6,5,7]))==gap_perm([2,3,1,4,6,5])
  @test_throws ArgumentError perm(G,[2,3,4,5,6,7,1])
  @test_throws ArgumentError perm(G, [1,1])
  @test one(G)==cperm(G,Int64[])
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

@testset "Eltypes" begin
   @test eltype(PermGroup)==PermGroupElem
   @test eltype(PcGroup)==PcGroupElem
   @test eltype(FPGroup)==FPGroupElem
   @test eltype(GL(2,3))==MatrixGroupElem{fq_nmod,fq_nmod_mat}
   @test eltype(DirectProductGroup)==Oscar.BasicGAPGroupElem{DirectProductGroup}
   @test eltype(direct_product(symmetric_group(3),cyclic_group(2)))==Oscar.BasicGAPGroupElem{DirectProductGroup}
   @test eltype(SemidirectProductGroup)==Oscar.BasicGAPGroupElem{SemidirectProductGroup}
   @test eltype(WreathProductGroup)==Oscar.BasicGAPGroupElem{WreathProductGroup}
   @test eltype(AutomorphismGroup{PcGroup})==Oscar.BasicGAPGroupElem{AutomorphismGroup{PcGroup}}

   G = symmetric_group(5)
   x = cperm([1,4,2,5])
   H = sub(G,[x])[1]
   y = cperm(G,[2,3,4])
   w = cperm(G,[1,4])
   K = sub(G,[w])[1]
   cc = conjugacy_class(G,y)
   cs = conjugacy_class(G,H)
   lc = x*H
   rc = H*x
   dc = K*x*H
   @test [z for z in G] == @inferred collect(G)
   @test [z for z in cc] == @inferred collect(cc)
   @test [z for z in cs] == @inferred collect(cs)
   @test [z for z in lc] == @inferred collect(lc)
   @test [z for z in rc] == @inferred collect(rc)
   @test [z for z in dc] == @inferred collect(dc)
   @test typeof(collect(G))==Vector{typeof(x)}
   @test typeof(collect(lc))==Vector{typeof(x)}
   @test typeof(collect(rc))==Vector{typeof(x)}
   @test typeof(collect(dc))==Vector{typeof(x)}
   @test typeof(collect(cc))==Vector{typeof(x)}
   @test typeof(collect(cs))==Vector{typeof(H)}
   @test eltype(cc)==typeof(y)
   @test eltype(cs)==typeof(H)
   @test eltype(lc)==typeof(y)
   @test eltype(rc)==typeof(y)
   @test eltype(dc)==typeof(y)
end

@testset "Generators" begin
   L=[symmetric_group(4), cyclic_group(5), free_group(3), symplectic_group(4,3)]
   for G in L
      K=gens(G)
      @test length(K) == ngens(G)
      for i in 1:length(K)
         @test K[i] == G[i]
         @test K[i] == gen(G,i)
      end
   end

   G = free_group(2)
   @test_throws AssertionError gen(G, 3)
   @test_throws ErrorException gen(G, 0)
end

@testset "deepcopy" begin
   for g in [symmetric_group(5), free_group(2), small_group(8, 1),
             automorphism_group(alternating_group(4))]
     m = Oscar.BasicGAPGroupElem(g, gen(g, 1).X)
     @test isdefined(m, :X)
     c = deepcopy(m);
     @test isdefined(c, :X)
     @test c.X == m.X

     @test deepcopy([one(g)]) == [one(g)]
   end
end
