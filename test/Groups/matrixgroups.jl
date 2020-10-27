
#FIXME : this may change in future. It can be easily skipped.
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

@testset "Iterator" begin
   G = SL(2,3)
   N = 0
   for x in G
      N+=order(x)
   end
   @test N==99
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
   @test_throws ArgumentError G([1,1,0,0])

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
   @test_throws ArgumentError O([z,0,0,1])

end

@testset "Methods on elements" begin
   T,t=PolynomialRing(FiniteField(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   x = G([1,z,0,1])
   y = G([z+2,0,0,1])
   @test x*y==G([z+2,z,0,1])
   @test x^-1==G([1,2*z,0,1])
   @test y^x==G([z+2,z+2,0,1])
   @test comm(y,x)==G([1,1,0,1])
   @test isone(G([1,0,0,1]))
   @test !isone(x)
   @test det(x)==1
   @test det(y)==z+2
   @test trace(x)==2
   @test tr(y)==z
   @test order(x)==3
   @test order(y)==8
   @test base_ring(x)==F
   @test nrows(y)==2
end

@testset "Subgroups" begin
   T,t=PolynomialRing(FiniteField(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   s1 = G([2,1,2,0])
   s2 = G([z+1,0,0,z+2])
   S,f = sub(G,[s1,s2])
   @test gens(S)==[s1,s2]
   @test base_ring(S)==F
   @test index(G,S)==8
   @test S==SL(2,F)
   @test parent(f(S[1]))==G
   @test f(S[1])==G(S[1])
   @test f(S[2])==G(S[2])
   O = GO(1,2,F)
   H = intersect(S,O)[1]
   @test H==SO(1,2,F)
   @test isnormal(O,H)
#   @test index(GO(0,3,3), omega_group(0,3,3))==4
#   @test index(GO(-1,4,2), omega_group(-1,4,2))==2
end

@testset "Cosets and conjugacy classes" begin
   T,t=PolynomialRing(FiniteField(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   H = GO(-1,2,F)
   x = G[1]
   lc = x*H
   @test order(lc)==order(H)
   @test representative(lc)==x
   @test acting_domain(lc)==H
   @test x in lc
   C = centralizer(G,x)[1]
   @test order(C)==64
   cc = conjugacy_class(G,x)
   @test x^G[2] in elements(cc)
   @test length(cc)==index(G,C)
end
