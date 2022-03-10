@testset "embeddings" begin
   @testset for G in [symmetric_group(5), small_group(24, 12), general_linear_group(2, 3)]
     G = symmetric_group(5)
     H, emb = sylow_subgroup(G, 2)
     x = gen(H, 1)
     y = image(emb, x)
     @test preimage(emb, y) == x
     @test any(g -> ! haspreimage(emb, g)[1], gens(G))
   end
end

n = 6
@testset "Homomorphism in Sym($n)" begin
   G = symmetric_group(n)
   x = G(vcat(2:n-1, [1]))

   g = hom(G, G, gens(G), [y^x for y in gens(G)])
   @test g == hom(G, G, gens(G), [y^x for y in gens(G)], check = false)
   @test g == hom(G, G, [y^x for y in gens(G)])
   @test g == hom(G, G, [y^x for y in gens(G)], check = false)
   @test typeof(g) == Oscar.GAPGroupHomomorphism{PermGroup, PermGroup}
   @test domain(g) == G
   @test codomain(g) == G
   @test image(g)[1] == G
   @test g(one(G)) == one(G)
   @test g(x)==x
   @test Set([y for y in centralizer(G,x)[1]]) == Set([y for y in G if g(y)==y])
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

   A = alternating_group(n)
   x = cperm(G,[1,2,3])
   f = hom(G,G,y -> y^x)
   @test typeof(restrict_homomorphism(f,A)) == Oscar.GAPGroupHomomorphism{PermGroup,PermGroup}
   fa = restrict_homomorphism(f,A)
   @test domain(fa)==A
   @test codomain(fa)==G
   @test fa(A[1])==f(A[1])
   @test_throws AssertionError restrict_homomorphism(fa, G)
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

   @testset "Change type" begin
       S = symmetric_group(4)
       (S1,f)=isomorphic_perm_group(S)
       @test S1==S
       (G,f)=isomorphic_pc_group(S)
       @test G isa PcGroup
       @test domain(f)==S
       @test codomain(f)==G
       @test isinjective(f)
       @test issurjective(f)

       (S,g)=isomorphic_pc_group(G)
       @test S isa PcGroup
       @test domain(g)==G
       @test codomain(g)==S
       @test isinjective(g)
       @test issurjective(g)

       (F,g)=isomorphic_fp_group(G)
       @test F isa FPGroup
       @test domain(g)==G
       @test codomain(g)==F
       @test isinjective(g)
       @test issurjective(g)

       @test_throws ArgumentError isomorphic_pc_group(symmetric_group(5))   
   end   
end

TestDirectProds=function(G1,G2)
   G = direct_product(G1,G2)
   f1 = embedding(G,1)
   f2 = embedding(G,2)
   p1 = projection(G,1)
   p2 = projection(G,2)

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
   G = direct_product(C2,C4)
   TestDirectProds(C2,C4)
   @test order(G)==8
   @test isabelian(G)
   @test !iscyclic(G)
   @test typeof(G)==DirectProductGroup
   @test Set([order(Int, x) for x in G])==Set([1,2,4])

   C3 = cyclic_group(3)
   C7 = cyclic_group(7)
   G = direct_product(C3,C7)
   TestDirectProds(C3,C7)
   @test order(G)==21
   @test isabelian(G)
   @test iscyclic(G)
   @test typeof(G)==DirectProductGroup
   @test Set([order(Int, x) for x in G])==Set([1,3,7,21])

   S4 = symmetric_group(4)
   A5 = alternating_group(5)
   G = direct_product(S4,A5)
   TestDirectProds(S4,A5)
   @test order(G)==1440
   @test typeof(G)==DirectProductGroup
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

@testset "Automorphism group of Sym(n)" begin
   G=symmetric_group(4)
   A=automorphism_group(G)

   @test A isa AutomorphismGroup
   @test A isa AutomorphismGroup{PermGroup}
   @test A.G == G
   @test isisomorphic(G,A)[1]
   @test order(A) == 24
   @test A==inner_automorphisms_group(A)[1]

   f = rand(A)
   g = rand(A)
   x = rand(G)
   o = order(f)
   fh = hom(f)
   @test f isa Oscar.GAPGroupElem{typeof(A)}
   @test fh isa Oscar.GAPGroupHomomorphism{PermGroup,PermGroup}
   @test A(fh)==f
   @test f(x)==x^f
   @test f^o == one(A)
   @test f*f^-1 == one(A)
   @test (f*g)(x) == g(f(x))
   @test comm(f,g) == f^-1*g^-1*f*g
   @test f(G[1])==fh(G[1])
   @test f(G[2])==fh(G[2])
   alt = alternating_group(4)
   N,e = sub(G,[alt[1],alt[2]])
   @test e*f==e*fh

   C=cyclic_group(2)
   g = hom(G,C,x -> C[1]^((1-sign(x))÷2) )
   @test f*g == fh*g
   @test kernel(f*g)==kernel(g)

   @test isinner_automorphism(f)
   g1 = inner_automorphism(G(alt[1]))
   @test !(g in A)
   g1 = A(g1)
   @test g1 in A
   g2 = A(inner_automorphism(G(alt[2])))
   AA,phi = sub(A,[g1,g2])
   @test isisomorphic(AA,alt)[1]
   @test index(A,AA)==2
   @test isnormal(A,AA)
   @test phi(AA[1])==AA[1]
   @test phi(AA[2])==AA[2]
   @test order(quo(A,AA)[1])==2
   @test isinvariant(f,alt)

   H = alternating_group(4)
   x = cperm(G,[1,2,3])
   f = A(hom(G,G,y->y^x))
   fa = restrict_automorphism(f,H)
   @test parent(fa)==automorphism_group(H)
   @testset for g in gens(H)
      @test fa(g)==H(f(g))
   end

   S = symmetric_group(3)
   g = hom(G,S,[cperm([1,2,3,4]), cperm([1,2])], [cperm([1,3]), cperm([1,2])])
   @test induced_automorphism(g,f)==automorphism_group(S)(inner_automorphism(cperm(S,[1,2,3])))
end

@testset "Other automorphisms groups" begin
   C = cyclic_group(3)
   G = direct_product(C,C)
   A = automorphism_group(G)

   @test isisomorphic(A,GL(2,3))[1]
   @test order(inner_automorphisms_group(A)[1])==1
end

@testset "Composition of mappings" begin
   g = symmetric_group(4)
   q, epi = quo(g, pcore(g, 2)[1])
   F, iso = isomorphic_perm_group(q)
   comp = compose(epi, iso)
   @test domain(comp) == domain(epi)
   @test codomain(comp) == codomain(iso)
end
