
n = 6
@testset "Homomorphism in Sym($n)" begin
   G = symmetric_group(n)
   x = rand(G)
   while x != one(G)
      x = rand(G)
   end

   g = hom(G, G, gens(G), [y^x for y in gens(G)])
   @test typeof(g) == Oscar.GAPGroupHomomorphism{PermGroup, PermGroup}
   @test domain(g) == G
   @test codomain(g) == G
   @test image(g)[1] == G
   @test g(one(G)) == one(G)
   @test g(x)==x
   @test Set(elements(centralizer(G,x)[1])) == Set([y for y in G if g(y)==y])
   z = rand(G)
   @test g(z) == z^x
   @test isinjective(g)
   @test issurjective(g)
   @test isinvertible(g)
   @test isbijective(g)
   og = order(g)
   @test og isa Integer
   @test og == order(x)
   @test g^og == id_hom(G)
   @test g^(og+1) == g
   @test g^(1-og) == g

   @test !isisomorphic(symmetric_group(4), symmetric_group(3))[1]
end

@testset "Operations on homomorphism in Sym($n)" begin
   G = symmetric_group(n)
   x = cperm(G,[1,2,3])
   Hx,fx = sub(G,[x])

   @test isinjective(fx)
   @test !issurjective(fx)
   @test fx(Hx[1])==x
   @test image(fx)[1] == Hx
   @test image(fx)[2] == fx

   y = cperm(G,[1,3,4])
   z = cperm(G,[1,2,5])
   Hy,fy = sub(G,[y])
   Hz,fz = sub(G,[z])
   f = hom(Hx,Hy,gens(Hx),gens(Hy))
   g = hom(Hy,Hz,gens(Hy),gens(Hz))
   @test isbijective(f)
   @test domain(f)==Hx
   @test codomain(f)==Hy
   @test image(f)[1]==Hy
   @test image(f)[2]==id_hom(Hy)
   @test f(x)==y
   @test inv(f)==f^-1
   @test domain(f^-1)==codomain(f)
   @test codomain(f^-1)==domain(f)
   @test image(f^-1)[1]==Hx
   @test image(f^-1)[2]==id_hom(Hx)
   @test (f^-1)(y)==x
   @test g(f(x))==z
   @test (f*g)(x)==z
   @test (f*g)^-1 == g^-1*f^-1
   @test_throws AssertionError g*f
   ty = trivial_morphism(Hy,Hy)
   @test f*ty==trivial_morphism(Hx,Hy)
   @test ty*g==trivial_morphism(Hy,Hz)
   @test f*trivial_morphism(Hy,Hz) == trivial_morphism(Hx,Hz)
   @test trivial_morphism(Hz,Hy)*f^-1 == trivial_morphism(Hz,Hx)
   @test order(kernel(f)[1])==1
   @test kernel(ty)[1] == Hy
   @test kernel(ty)[2] == id_hom(Hy)
   @test kernel(f*ty)[1]==Hx
   
end

@testset "Isomorphic groups" begin
   @testset "Dihedral_as_permutation" for n in 4:10
      G = symmetric_group(n)
      D = dihedral_group(PermGroup,2*n)
      H = sub(G, [G(D[1]),G(D[2])])[1]
   
      @test order(H)==2*n
      @test H == D
   end

   @testset "Abelian_as_permutation" for n in 15:20
      G = symmetric_group(n)
      for j in 2:n-2
         D = abelian_group(PermGroup,[j,n-j])
         H = sub(G, [G(D[1]),G(D[2])])[1]

         @test H == D
         @test order(D) == j*(n-j)
         @test (gcd(j,n-j)==1 && iscyclic(D)) || (gcd(j,n-j)>1 && !iscyclic(D))
      end
   end

end

TestDirectProds=function(G1,G2)
   G,f1,f2,p1,p2 = direct_product(G1,G2, :both)

   @test isinjective(f1)
   @test isinjective(f2)
   @test issurjective(p1)
   @test issurjective(p2)
   @test domain(f1)==G1
   @test codomain(f1)==G
   @test domain(f2)==G2
   @test codomain(f2)==G
   @test domain(p1)==G
   @test codomain(p1)==G1
   @test domain(p2)==G
   @test codomain(p2)==G2
   @test f1*p1==id_hom(G1)
   @test f2*p2==id_hom(G2)
   for i in 1:ngens(G1)
      @test f1(G1[i])==G[i]
   end
   for i in 1:ngens(G2)
      @test f2(G2[i])==G[i+ngens(G1)]
   end
   q1=p1*f1
   q2=p2*f2
   @test isisomorphic(kernel(q1)[1],G2)[1]
   @test isisomorphic(image(q1)[1],G1)[1]
   @test isisomorphic(kernel(q2)[1],G1)[1]
   @test isisomorphic(image(q2)[1],G2)[1]
end

@testset "Direct product" begin
   C2 = cyclic_group(2)
   C4 = cyclic_group(4)
   G = direct_product(C2,C4)[1]
   TestDirectProds(C2,C4)
   @test order(G)==8
   @test isabelian(G)
   @test !iscyclic(G)
   @test typeof(G)==PcGroup
   @test Set([order(x) for x in G])==Set([1,2,4])

   C3 = cyclic_group(3)
   C7 = cyclic_group(7)
   G = direct_product(C3,C7)[1]
   TestDirectProds(C3,C7)
   @test order(G)==21
   @test isabelian(G)
   @test iscyclic(G)
   @test typeof(G)==PcGroup
   @test Set([order(x) for x in G])==Set([1,3,7,21])

   S4 = symmetric_group(4)
   A5 = alternating_group(5)
   G = direct_product(S4,A5)[1]
   TestDirectProds(S4,A5)
   @test order(G)==1440
   @test typeof(G)==PermGroup
end

TestKernels = function(G,H,f)
   K,i = kernel(f)
   Im = image(f)[1]

   @test preimage(f,H)==(G,id_hom(G))
   @test preimage(f,sub(H,[one(H)])[1])==(K,i)
   z=rand(Im)
   @test haspreimage(f,z)[1]
   @test f(haspreimage(f,z)[2])==z

   @test isinjective(i)
   for j in 1:ngens(K)
      @test i(K[j]) in G
      @test (i*f)(K[j])==one(H)
      @test index(G,K)==order(Im)
   end
   if isnormal(H,Im) 
      C,p = cokernel(f)
        @test issurjective(p)
      for j in 1:ngens(G)
         @test (f*p)(G[j])==one(C)
      end
      @test index(H,Im)==order(C)
   end
end
   
@testset "Kernel and cokernel" begin
   G=symmetric_group(4)
   z=rand(G)
   TestKernels(G,G,hom(G,G, x -> x^z))
   
   C=cyclic_group(2)
   TestKernels(G,C,hom(G,C,gens(G),[C[1],C[1]]))      #sign

   G=GL(2,7)
   C=cyclic_group(6)
   TestKernels(G,C,hom(G,C,gens(G),[C[1],one(C)]))        #determinant

   G=abelian_group(PcGroup,[3,3,3])
   H=abelian_group(PcGroup,[3,3])
   f=hom(G,H,gens(G),[H[1],one(H),one(H)])
   TestKernels(G,H,f)
end

