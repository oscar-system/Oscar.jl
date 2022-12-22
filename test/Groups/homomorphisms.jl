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
   @test is_injective(g)
   @test is_surjective(g)
   @test is_invertible(g)
   @test is_bijective(g)
   og = order(g)
   @test og isa Integer
   @test og == order(x)
   @test g^og == id_hom(G)
   @test g^(og+1) == g
   @test g^(1-og) == g

   @test !is_isomorphic(symmetric_group(4), symmetric_group(3))

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

   @test is_injective(fx)
   @test !is_surjective(fx)
   @test fx(Hx[1])==x
   @test image(fx)[1] == Hx
   @test image(fx)[2] == fx

   y = cperm(G,[1,3,4])
   z = cperm(G,[1,2,5])
   Hy,fy = sub(G,[y])
   Hz,fz = sub(G,[z])
   f = hom(Hx,Hy,gens(Hx),gens(Hy))
   g = hom(Hy,Hz,gens(Hy),gens(Hz))
   @test is_bijective(f)
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

@testset "quo for FPGroup" begin
   # quotient of a *full* free group
   g = free_group(2)
   x, y = gens(g)
   grels = [x^2, y^2, comm(x, y)]
   q, epi = quo(g, grels)
   @test ngens(q) == 2
   @test order(q) == 4
   @test [epi(h) for h in gens(g)] == gens(q)

   # quotient of a *full* f.p. group
   q2, epi2 = quo(q, [q[1]])
   @test order(q2) == 2
   @test [epi2(h) for h in gens(q)] == gens(q2)

   # quotient of a *subgroup* of a free group:
   # We forbid this on the Oscar side,
   # otherwise GAP would first construct the quotient group
   # but later run into an error when one asks for its order.
   h = sub(g, [x^2])[1]
   @test_throws AssertionError quo(h, [x^10])
   # q3, epi3 = quo(h, [x^10])
   # @test order(q3) == 5
   # @test [epi3(h) for h in gens(h)] == gens(q3)

   # quotient of a *subgroup* of a f.p. group:
   # We forbid this on the Oscar side,
   # otherwise GAP would run into an error.
   h = sub(q2, [q2[2]^5])[1]
   @test_throws AssertionError quo(h, [h[1]^2])
end

@testset "map_word" begin
   # Create a free group in GAP in syllable words family,
   # in order to make the tests.
   GAP.Globals.PushOptions(GAP.GapObj(Dict(:FreeGroupFamilyType => GAP.GapObj("syllable"))))
   FS = free_group(2)   # syllable representation
   GAP.Globals.PopOptions()
   FL = free_group(2)   # letter representation
   @test GAP.Globals.IsSyllableWordsFamily(
           GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(FS.X)))
   @test GAP.Globals.IsLetterWordsFamily(
           GAP.Globals.ElementsFamily(GAP.Globals.FamilyObj(FL.X)))
   for F in [FL, FS]
     F1, F2 = gens(F)
     rels = [F1^2, F2^2, comm(F1, F2)]
     FP = quo(F, rels)[1]

     # map an element of a free group or a f.p. group ...
     for f in [F, FP]
       f1, f2 = gens(f)
       of = one(f)

       # ... to an element of a permutation group or to a rational number
       for imgs in [gens(symmetric_group(4)), fmpq[2, 3]]
         g1 = imgs[1]
         g2 = imgs[2]
         for (x, w) in [(f1, g1), (f2, g2), (f2^2*f1^-3, g2^2*g1^-3)]
           @test map_word(x, imgs) == w
         end
         @test map_word(of, imgs, init = 0) == one(g1)  # `init` is ignored
       end
     end

     # empty list of generators
     T = free_group(0)
     @test map_word(one(T), [], init = 1) == 1            # `init` is returned
     @test map_word(one(F), gens(F), init = 1) == one(F)  # `init` is ignored

     # wrong number of images
     @test_throws AssertionError map_word(one(F), [])
     @test_throws AssertionError map_word(one(F), [F1])
     @test_throws AssertionError map_word(one(F), [F1, F1, F1])
   end

   # map according to a description of the word
   for imgs in [gens(symmetric_group(4)), fmpq[2, 3]]
     g1 = imgs[1]
     g2 = imgs[2]
     for (v, w) in [
       # via exponents
       ([-1, 2, -1, 2], (g1^-1*g2)^2),
       ([-1, -1, -1, 2, 2, 2], g1^-3*g2^3),
       # via pairs
       ([1 => -1, 2 => 1, 1 => -1, 2 => 1], (g1^-1*g2)^2),
       ([1 => -3, 2 => 3], g1^-3*g2^3),
       ]

       @test map_word(v, imgs) == w
       invs = Vector(undef, 2)
       @test map_word(v, imgs, genimgs_inv = invs) == w
       @test isassigned(invs, 1)
       @test ! isassigned(invs, 2)
     end
   end

   # empty list of generators
   @test map_word([], [], init = 0) == 0        # `init` is returned
   @test map_word([], [2, 3], init = 0) == 1    # `init` is ignored

   # wrong number of images
   @test_throws AssertionError map_word([3], [])
   @test_throws AssertionError map_word([-3], [2, 3])
   @test_throws AssertionError map_word([3 => 1], [2, 3])
end

@testset "Isomorphic groups" begin
   @testset "Dihedral_as_permutation" for n in 4:10
      G = symmetric_group(n)
      D = dihedral_group(PermGroup,2*n)
      H = sub(G, [G(D[1]),G(D[2])])[1]
   
      @test order(H)==2*n
      @test H == D
   end

   @testset "Finite abelian GAPGroup to GrpAbFinGen" begin
      for invs in [[1], [2, 3, 4], [6, 8, 9, 15]], T in [PermGroup, PcGroup, FPGroup]
         G = abelian_group(T, invs)
         iso = @inferred isomorphism(GrpAbFinGen, G)
         A = codomain(iso)
         @test order(G) == order(A)
         for x in gens(G)
            y = image(iso, x)
            @test preimage(iso, y) == x
         end
      end
   end

   @testset "Finite GrpAbFinGen to GAPGroup" begin
      @testset for Agens in [[2, 4, 8], [2, 3, 4], [ 2, 12 ],
                             [1, 6], matrix(ZZ, 2, 2, [2, 3, 2, 6])]
         A = abelian_group(Agens)
         for T in [FPGroup, PcGroup, PermGroup]
            iso = @inferred isomorphism(T, A)
            for x in gens(A), y in gens(A)
               z = x+y
               @test iso(x) * iso(y) == iso(z)
               @test all(a -> preimage(iso, iso(a)) == a, [x, y, z])
            end
         end
      end
   end

   @testset "Infinite GrpAbFinGen to GAPGroup" begin
      Agens = matrix(ZZ, 2, 2, [2, 3, 0, 0])
      A = abelian_group(Agens)
      for T in [FPGroup]
         iso = @inferred isomorphism(T, A)
         for x in gens(A), y in gens(A)
            z = x+y
            @test iso(x) * iso(y) == iso(z)
            @test all(a -> preimage(iso, iso(a)) == a, [x, y, z])
         end
      end
   end

   @testset "GrpAbFinGen to GrpAbFinGen" begin
      A = abelian_group([2, 3, 4])
      iso = @inferred isomorphism(GrpAbFinGen, A)
   end

   @testset "GrpGen to GAPGroups" begin
      G = Hecke.small_group(64, 14, DB = Hecke.DefaultSmallGroupDB())
      for T in [FPGroup, PcGroup, PermGroup]
         iso = @inferred isomorphism(T, G)
         for x in gens(G), y in gens(G)
            z = x * y
            @test iso(x) * iso(y) == iso(z)
            @test all(a -> preimage(iso, iso(a)) == a, [x, y, z])
         end
      end

      H = small_group(64, 14)
      @test isisomorphic(G, H)
      f = isomorphism(G, H)
      for x in gens(G), y in gens(G)
         @test f(x) * f(y) == f(x * y)
         @test preimage(f, f(x)) == x
         @test preimage(f, f(y)) == y
      end
      fl, f = is_isomorphic_with_map(G, H)
      @test fl
      for x in gens(G), y in gens(G)
         @test f(x) * f(y) == f(x * y)
         @test preimage(f, f(x)) == x
         @test preimage(f, f(y)) == y
      end

      @test isisomorphic(H, G)
      f = isomorphism(H, G)
      for x in gens(H), y in gens(H)
         @test f(x) * f(y) == f(x * y)
         @test preimage(f, f(x)) == x
         @test preimage(f, f(y)) == y
      end
      fl, f = is_isomorphic_with_map(H, G)
      @test fl
      for x in gens(H), y in gens(H)
         @test f(x) * f(y) == f(x * y)
         @test preimage(f, f(x)) == x
         @test preimage(f, f(y)) == y
      end

      H = cyclic_group(2)
      @test !isisomorphic(G, H)
      @test_throws ArgumentError isomorphism(G, H)
      fl, _ = is_isomorphic_with_map(G, H)
      @test !fl
      @test !isisomorphic(H, G)
      @test_throws ArgumentError isomorphism(H, G)
      fl, _ = is_isomorphic_with_map(H, G)
      @test !fl
   end

   @testset "Group types as constructors" begin
      G = symmetric_group(4)
      for T in [FPGroup, PcGroup, PermGroup]
        H = T(G)
        @test H isa T
        @test is_isomorphic(G, H)[1]
      end

      G = cyclic_group(5)
      T = GrpAbFinGen
      H = T(G)
      @test H isa T
      @test order(H) == order(G)
      K = PermGroup(H)
      @test K isa PermGroup
      @test order(K) == order(H)
   end

   @testset "Abelian_as_permutation" for n in 15:20
      G = symmetric_group(n)
      for j in 2:n-2
         D = abelian_group(PermGroup,[j,n-j])
         H = sub(G, [G(D[1]),G(D[2])])[1]

         @test H == D
         @test order(D) == j*(n-j)
         @test (gcd(j,n-j)==1 && is_cyclic(D)) || (gcd(j,n-j)>1 && !is_cyclic(D))
      end
   end

   @testset "Change type" begin
       S = symmetric_group(4)
       f = @inferred isomorphism(PermGroup, S)
       @test codomain(f) == S

       f = @inferred isomorphism(PcGroup, S)
       G = codomain(f)
       @test G isa PcGroup
       @test domain(f) == S
       @test is_injective(f)
       @test is_surjective(f)

       f = @inferred isomorphism(PcGroup, G)
       @test codomain(f) isa PcGroup
       @test domain(f) == G
       @test is_injective(f)
       @test is_surjective(f)

       f = @inferred isomorphism(FPGroup, G)
       @test codomain(f) isa FPGroup
       @test domain(f) == G
       @test is_injective(f)
       @test is_surjective(f)

       G = abelian_group(PermGroup, [2, 2])
       f = @inferred isomorphism(GrpAbFinGen, G)
       @test codomain(f) isa GrpAbFinGen
       @test domain(f) == G
     # @test is_injective(f)
     # @test is_surjective(f)

       @test_throws ArgumentError isomorphism(GrpAbFinGen, symmetric_group(5))
       @test_throws ArgumentError isomorphism(PcGroup, symmetric_group(5))
       @test_throws ArgumentError isomorphism(PermGroup, free_group(1))

       G = symmetric_group(4)
       @test PermGroup(G) isa PermGroup
       @test PcGroup(G) isa PcGroup
       @test FPGroup(G) isa FPGroup
       @test_throws ArgumentError GrpAbFinGen(G)
   end
end

TestDirectProds=function(G1,G2)
   G = direct_product(G1,G2)
   f1 = embedding(G,1)
   f2 = embedding(G,2)
   p1 = projection(G,1)
   p2 = projection(G,2)

   @test is_injective(f1)
   @test is_injective(f2)
   @test is_surjective(p1)
   @test is_surjective(p2)
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
   @test is_isomorphic(kernel(q1)[1],G2)
   @test is_isomorphic(image(q1)[1],G1)
   @test is_isomorphic(kernel(q2)[1],G1)
   @test is_isomorphic(image(q2)[1],G2)
end

@testset "Direct product" begin
   C2 = cyclic_group(2)
   C4 = cyclic_group(4)
   G = direct_product(C2,C4)
   TestDirectProds(C2,C4)
   @test order(G)==8
   @test is_abelian(G)
   @test !is_cyclic(G)
   @test typeof(G)==DirectProductGroup
   @test Set([order(Int, x) for x in G])==Set([1,2,4])

   C3 = cyclic_group(3)
   C7 = cyclic_group(7)
   G = direct_product(C3,C7)
   TestDirectProds(C3,C7)
   @test order(G)==21
   @test is_abelian(G)
   @test is_cyclic(G)
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

   @test is_injective(i)
   for j in 1:ngens(K)
      @test i(K[j]) in G
      @test (i*f)(K[j])==one(H)
      @test index(G,K)==order(Im)
   end
   if is_normal(H,Im) 
      C,p = cokernel(f)
        @test is_surjective(p)
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
   @test is_isomorphic(G,A)
   @test order(A) == 24
   @test A==inner_automorphism_group(A)[1]

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

   @test is_inner_automorphism(f)
   g1 = inner_automorphism(G(alt[1]))
   @test !(g in A)
   g1 = A(g1)
   @test g1 in A
   g2 = A(inner_automorphism(G(alt[2])))
   AA,phi = sub(A,[g1,g2])
   @test is_isomorphic(AA,alt)
   @test index(A,AA)==2
   @test is_normal(A,AA)
   @test phi(AA[1])==AA[1]
   @test phi(AA[2])==AA[2]
   @test order(quo(A,AA)[1])==2
   @test is_invariant(f,alt)

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

   @test is_isomorphic(A,GL(2,3))
   @test order(inner_automorphism_group(A)[1])==1
end

@testset "Composition of mappings" begin
   g = symmetric_group(4)
   q, epi = quo(g, pcore(g, 2)[1])
   iso = @inferred isomorphism(PermGroup, q)
   comp = compose(epi, iso)
   @test domain(comp) == domain(epi)
   @test codomain(comp) == codomain(iso)
end
