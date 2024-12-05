@testset "embeddings" begin
   @testset for G in [symmetric_group(5), small_group(24, 12), general_linear_group(2, 3)]
     G = symmetric_group(5)
     H, emb = sylow_subgroup(G, 2)
     x = gen(H, 1)
     y = image(emb, x)
     @test preimage(emb, y) == x
     @test any(g -> ! has_preimage_with_preimage(emb, g)[1], gens(G))
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
   # quotient of a *full* free group by a vector of elements
   g = free_group(2)
   x, y = gens(g)
   grels = [x^2, y^2, comm(x, y)]
   q, epi = quo(g, grels)
   @test ngens(q) == 2
   @test order(q) == 4
   @test [epi(h) for h in gens(g)] == gens(q)

   # quotient of a *full* f.p. group by a vector of elements
   # This is handled via GAP's `\/`.
   q2, epi2 = quo(q, [q[1]])
   @test order(q2) == 2
   @test [epi2(h) for h in gens(q)] == gens(q2)

   # quotient of a *subgroup* of a free or f.p. group by a vector of elements
   # We forbid this call on the Oscar side.
   # (Note:
   # - Currently GAP's `G / elms` cannot handle this;
   #   up to GAP 4.12.2, GAP first constructs a quotient group
   #   but later runs into an error when one asks for its order.
   # - Trying to compute a quotient by the normal closure
   #   in `G` of the given `elms` may or may not run into a GAP error.)
   h = sub(g, [x^2])[1]
   @test_throws ArgumentError quo(h, [h(x^10)])
   n = normal_closure(h, sub(h, [h(x^10)])[1])[1]
   @test_throws ErrorException quo(h, n)

   h = sub(q2, [q2[2]^5])[1]
   @test_throws ArgumentError quo(h, [h[1]^2])
end

@testset "map_word for f.p. groups" begin
   # Create a free group in GAP in syllable words family,
   # in order to make the tests.
   GAP.Globals.PushOptions(GapObj(Dict(:FreeGroupFamilyType => GapObj("syllable"))))
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

     # map an element of a (subgroup of a) free group or a f.p. group ...
     for f in [F, FP]
       f1, f2 = gens(f)
       of = one(f)
       s1 = gen(sub(f, [f1])[1], 1)

       # ... to an element of a permutation group or to a rational number ...
       for imgs in [gens(symmetric_group(4)),
                    QQFieldElem[2, 3]]
         g1 = imgs[1]
         g2 = imgs[2]
         for (x, w) in [(f1, g1), (f2, g2), (f2^2*f1^-3, g2^2*g1^-3),
                        (s1^2, g1^2)]
           @test map_word(x, imgs) == w
           @test map_word(x, imgs, init = g1) == g1*w  # `init` is used
         end
         @test map_word(of, imgs, init = 0) == 0  # `init` is returned
       end

       # ... or to an element in an (additive) abelian group
       imgs = gens(abelian_group(3, 5))
       g1 = imgs[1]
       g2 = imgs[2]
       for (x, w) in [(f1, g1), (f2, g2), (f2^2*f1^-3, 2*g2-3*g1)]
         @test map_word(x, imgs) == w
         @test map_word(x, imgs, init = g1) == g1+w  # `init` is used
       end
       @test map_word(of, imgs, init = 0) == 0  # `init` is returned
     end

     # empty list of generators
     T = free_group(0)
     @test_throws ArgumentError map_word(one(T), Int[])  # no `init`
     @test map_word(one(T), [], init = 1) == 1           # `init` is returned
     @test map_word(one(F), gens(F)) == one(F)           # works without `init`

     # wrong number of images
     @test_throws AssertionError map_word(one(F), [])
     @test_throws AssertionError map_word(one(F), [F1])
     @test_throws AssertionError map_word(one(F), [F1, F1, F1])
   end

   # map according to a description of the word
   for imgs in [gens(symmetric_group(4)), QQFieldElem[2, 3]]
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
   @test map_word([], [2, 3], init = 0) == 0    # `init` is returned
   @test map_word([], [2, 3]) == 1              # no `init` given, try `one`
   G = abelian_group(2)
   @test map_word([], gens(G)) == zero(G)

   # wrong number of images
   @test_throws AssertionError map_word([3], [])
   @test_throws AssertionError map_word([-3], [2, 3])
   @test_throws AssertionError map_word([3 => 1], [2, 3])
end

@testset "map_word for (sub) pc groups" begin
   for G in [ PcGroup(symmetric_group(4)),         # GAP Pc Group
            # abelian_group(PcGroup, [2, 3, 4]),   # problem with gens vs. pcgs
              abelian_group(PcGroup, [0, 3, 4]) ]  # GAP Pcp group
     n = number_of_generators(G)
     F = free_group(n)
     S = sub(G, gens(G))[1]
     for g in [G, S]
       for x in [one(g), rand(g)]
         img = map_word(x, gens(F))
         @test x == map_word(img, gens(g))
         invs = Vector(undef, n)
         img = map_word(x, gens(F), genimgs_inv = invs)
         @test x == map_word(img, gens(g))
       end
     end
   end
end

@testset "Isomorphic groups" begin
   @testset "Dihedral_as_permutation" for n in 4:10
      G = symmetric_group(n)
      D = dihedral_group(PermGroup,2*n)
      H = sub(G, [G(D[1]),G(D[2])])[1]
   
      @test order(H)==2*n
      @test H == D
   end

   @testset "Finite abelian GAPGroup to FinGenAbGroup" begin
#     for invs in [[1], [2, 3, 4], [6, 8, 9, 15]], T in [PermGroup, PcGroup, FPGroup]
      for invs in [[1], [2, 3, 4], [6, 8, 9, 15]], T in [PermGroup, SubPcGroup, FPGroup]
         G = abelian_group(T, invs)
         iso = @inferred isomorphism(FinGenAbGroup, G)
         A = codomain(iso)
         @test order(G) == order(A)
         for x in gens(G)
            y = image(iso, x)
            @test preimage(iso, y) == x
         end
      end
   end

   @testset "Finite FinGenAbGroup to GAPGroup" begin
      @testset for Agens in [Int[], [2, 4, 8], [2, 3, 4], [2, 12],
                             [1, 6], matrix(ZZ, 2, 2, [2, 3, 2, 6])]
         A = abelian_group(Agens)
         for T in [FPGroup, PcGroup, SubPcGroup, PermGroup]
            iso = @inferred isomorphism(T, A)
            for x in gens(A), y in gens(A)
               z = x+y
               @test iso(x) * iso(y) == iso(z)
               @test all(a -> preimage(iso, iso(a)) == a, [x, y, z])
            end
         end
      end
   end

   @testset "Infinite FinGenAbGroup to GAPGroup" begin
      @testset for Agens in [matrix(ZZ, 2, 2, [2, 3, 0, 0]), [6, 0]]
         A = abelian_group(Agens)
         for T in [FPGroup, PcGroup]
            iso = @inferred isomorphism(T, A)
            for x in gens(A), y in gens(A)
               z = x+y
               @test iso(x) * iso(y) == iso(z)
               @test all(a -> preimage(iso, iso(a)) == a, [x, y, z])
            end
         end
      end
   end

   @testset "FinGenAbGroup to FinGenAbGroup" begin
      A = abelian_group([2, 3, 4])
      iso = @inferred isomorphism(FinGenAbGroup, A)
   end

   @testset "Vector space to FPGroup or PcGroup" begin
      for (F, n) in [(GF(2), 0), (GF(3), 2), (GF(4), 2),
                     (Nemo.Native.GF(2), 0),
                     (Nemo.Native.GF(3), 2)]
        V = free_module(F, n)
        for T in [FPGroup, PcGroup]
          iso = @inferred isomorphism(T, V)
          for x in [zero(V), rand(V)]
            @test preimage(iso, iso(x)) == x
          end
        end
      end
   end

   @testset "MultTableGroup to GAPGroups" begin
      for G in [Hecke.small_group(64, 14, DB = Hecke.DefaultSmallGroupDB()),
                Hecke.small_group(20, 3, DB = Hecke.DefaultSmallGroupDB())]
         for T in [FPGroup, PcGroup, PermGroup]
            iso = @inferred isomorphism(T, G)
            for x in gens(G), y in gens(G)
               z = x * y
               @test iso(x) * iso(y) == iso(z)
               @test all(a -> preimage(iso, iso(a)) == a, [x, y, z])
            end
         end
      end   

      G  = Hecke.small_group(64, 14, DB = Hecke.DefaultSmallGroupDB())
      H = small_group(64, 14)
      @test is_isomorphic(G, H)
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

      @test is_isomorphic(H, G)
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
      @test !is_isomorphic(G, H)
      @test_throws ArgumentError isomorphism(G, H)
      fl, _ = is_isomorphic_with_map(G, H)
      @test !fl
      @test !is_isomorphic(H, G)
      @test_throws ArgumentError isomorphism(H, G)
      fl, _ = is_isomorphic_with_map(H, G)
      @test !fl
   end

   @testset "Group types as constructors" begin
      @testset "Source $G" for G in [
            cyclic_group(5),
            dihedral_group(10),
            symmetric_group(4),
            transitive_group(5,2),
            #abelian_group(5),  # FIXME error in is_isomorphic
            ]
         @testset "Range type $T" for (T, f) in [
              (FPGroup, fp_group),
              (PcGroup, pc_group),
              (PermGroup, permutation_group),
              #(FinGenAbGroup, FinGenAbGroup),  # FIXME: errors
              ]
            H = T(G)
            @test H isa T
            @test has_order(H)
            @test is_isomorphic(G, H)[1]

            H = f(G)
            @test H isa T
            @test has_order(H)
            @test is_isomorphic(G, H)[1]
         end
      end

      G = cyclic_group(5)
      T = FinGenAbGroup
      H = T(G)
      @test H isa T
      @test order(H) == order(G)
      K = PermGroup(H)
      @test K isa PermGroup
      @test order(K) == order(H)
      K = permutation_group(H)
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

       @test_throws ArgumentError isomorphism(PcGroup, S, on_gens = true)
       f = isomorphism(PcGroup, G, on_gens = true)
       @test [f(x) for x in gens(G)] == gens(codomain(f))

       f = @inferred isomorphism(PcGroup, G)
       @test codomain(f) isa PcGroup
       @test domain(f) == G
       @test is_injective(f)
       @test is_surjective(f)

       G = symmetric_group(5)
       f = @inferred isomorphism(FPGroup, G)
       @test codomain(f) isa FPGroup
       @test domain(f) == G
       @test is_injective(f)
       @test is_surjective(f)

       f2 = @inferred isomorphism(FPGroup, G, on_gens=true)
       @test codomain(f2) isa FPGroup
       @test domain(f2) == G
       @test is_injective(f2)
       @test is_surjective(f2)
       @test [preimage(f2, x) for x in gens(codomain(f2))] == gens(G)
       @test [preimage(f, x) for x in gens(codomain(f))] != gens(G)

       G = symmetric_group(1)
       iso = @inferred isomorphism(FPGroup, G, on_gens = true)
       @test ngens(G) == ngens(codomain(iso))
       @test is_bijective(iso)
       G = sub(G, [one(G)])[1]
       iso = @inferred isomorphism(FPGroup, G, on_gens = true)
       @test ngens(G) == ngens(codomain(iso))

       G = abelian_group(PermGroup, [2, 2])
       f = @inferred isomorphism(FinGenAbGroup, G)
       @test codomain(f) isa FinGenAbGroup
       @test domain(f) == G
     # @test is_injective(f)   # no method for GroupIsomorphismFromFunc
     # @test is_surjective(f)

       @test_throws ArgumentError isomorphism(FinGenAbGroup, symmetric_group(5))
       @test_throws ArgumentError isomorphism(PcGroup, symmetric_group(5))
       @test_throws ArgumentError isomorphism(PermGroup, free_group(1))

       G = symmetric_group(4)
       @test PermGroup(G) isa PermGroup
       @test permutation_group(G) isa PermGroup
       @test pc_group(G) isa PcGroup
       @test FPGroup(G) isa FPGroup
       @test_throws ArgumentError FinGenAbGroup(G)
   end
end

@testset "Homomorphism GAPGroup to FinGenAbGroup" begin
   # G abelian, A isomorphic to G
   G = abelian_group( PermGroup, [ 2, 4 ] )
   A = abelian_group( [ 2, 4 ] )
   imgs = gens(A)
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 1

   # G abelian, A a proper factor of G
   G = abelian_group( PermGroup, [ 2, 4 ] )
   A = abelian_group( [ 2, 2 ] )
   imgs = gens(A)
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 2

   # G abelian, A containing a proper factor of G
   G = abelian_group( PermGroup, [ 2, 4 ] )
   A = abelian_group( [ 2, 4 ] )
   imgs = [gen(A, 1), 2*gen(A, 2)]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 2

   # G nonabelian, A isomorphic to G/G'
   G = dihedral_group(8)
   A = abelian_group( [ 2, 2 ] )
   imgs = [gen(A, 1), gen(A, 2), zero(A)]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 2

   # G nonabelian, A a proper factor of  G/G'
   G = dihedral_group(8)
   A = abelian_group( [ 2 ] )
   imgs = [gen(A, 1), gen(A, 1), zero(A)]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 4

   # G nonabelian, A containing a proper factor of G/G'
   G = dihedral_group(8)
   A = abelian_group( [ 4 ] )
   imgs = [2*gen(A,1), 2*gen(A,1), zero(A)]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 4

   # G trivial
   G = cyclic_group(PcGroup, 1)
   A = abelian_group( [ 2 ] )
   imgs = elem_type(A)[]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 1

   # A trivial
   G = dihedral_group(8)
   A = abelian_group( [ 1 ] )
   imgs = [zero(A), zero(A), zero(A)]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 8

   # G and A trivial
   G = cyclic_group(PcGroup, 1)
   A = abelian_group( [ 1 ] )
   imgs = elem_type(A)[]
   mp = hom(G, A, imgs)
   @test order(kernel(mp)[1]) == 1
end

function test_direct_prods(G1,G2)
   G = direct_product(G1,G2)
   f1 = canonical_injection(G,1)
   f2 = canonical_injection(G,2)
   p1 = canonical_projection(G,1)
   p2 = canonical_projection(G,2)

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
   test_direct_prods(C2,C4)
   @test order(G)==8
   @test is_abelian(G)
   @test !is_cyclic(G)
   @test typeof(G)==DirectProductGroup
   @test Set([order(Int, x) for x in G])==Set([1,2,4])

   C3 = cyclic_group(3)
   C7 = cyclic_group(7)
   G = direct_product(C3,C7)
   test_direct_prods(C3,C7)
   @test order(G)==21
   @test is_abelian(G)
   @test is_cyclic(G)
   @test typeof(G)==DirectProductGroup
   @test Set([order(Int, x) for x in G])==Set([1,3,7,21])

   S4 = symmetric_group(4)
   A5 = alternating_group(5)
   G = direct_product(S4,A5)
   test_direct_prods(S4,A5)
   @test order(G)==1440
   @test typeof(G)==DirectProductGroup
end

function test_kernel(G,H,f)
   K,i = kernel(f)
   Im = image(f)[1]

#TODO: activate these tests as soon as they pass again;
#      the point is that comparing the embeddings is done via `===`
#  @test preimage(f,H)==(G,id_hom(G))
   @test preimage(f,sub(H,[one(H)])[1])==(K,i)
   z=rand(Im)
   @test has_preimage_with_preimage(f,z)[1]
   @test f(has_preimage_with_preimage(f,z)[2])==z

   @test is_injective(i)
   for j in 1:ngens(K)
      @test i(K[j]) in G
      @test (i*f)(K[j])==one(H)
      @test index(G,K)==order(Im)
   end
   if is_normalized_by(Im, H)
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
   test_kernel(G,G,hom(G,G, x -> x^z))
   
   C=cyclic_group(2)
   test_kernel(G,C,hom(G,C,gens(G),[C[1],C[1]]))      #sign

   G=GL(2,7)
   C=cyclic_group(6)
   test_kernel(G,C,hom(G,C,gens(G),[C[1],one(C)]))        #determinant

   G=abelian_group(PcGroup,[3,3,3])
   H=abelian_group(PcGroup,[3,3])
   f=hom(G,H,gens(G),[H[1],one(H),one(H)])
   test_kernel(G,H,f)
end

@testset "Automorphism group of a perm. group or a (sub) pc group" begin
   for T in [PermGroup, PcGroup, SubPcGroup]
      G = small_group(T, 24, 12)
      A = automorphism_group(G)

      @test A isa AutomorphismGroup
      @test A isa AutomorphismGroup{T}
      @test A.G === G
      @test is_isomorphic(G, A)
      @test order(A) == 24
      @test A == inner_automorphism_group(A)[1]

      f = rand(A)
      g = rand(A)
      x = rand(G)
      o = order(f)
      fh = hom(f)
      @test f isa Oscar.GAPGroupElem{typeof(A)}
      @test fh isa Oscar.GAPGroupHomomorphism{T, T}
      @test A(fh) == f
      @test f(x) == x^f
      @test f^o == one(A)
      @test f*f^-1 == one(A)
      @test (f*g)(x) == g(f(x))
      @test comm(f, g) == f^-1*g^-1*f*g
      @test f(G[1]) == fh(G[1])
      @test f(G[2]) == fh(G[2])
      H = derived_subgroup(G)[1]
      N, e = sub(G, [H[1], H[2]])
      @test e*f == e*fh

      C = cyclic_group(2) # type is independent of `T`
      oC = one(C)
      gC = gen(C, 1)
      g = hom(G, C, x -> x in H ? oC : gC)
      @test f*g == fh*g
      @test kernel(f*g) == kernel(g)
      @test induced_automorphism(g, f) == induced_automorphism(g, fh)

      @test is_inner_automorphism(f)
      g1 = inner_automorphism(G(H[1]))
      @test !(g in A)
      g1 = A(g1)
      @test g1 in A
      g2 = A(inner_automorphism(G(H[2])))
      AA, phi = sub(A, [g1, g2])
      @test is_isomorphic(AA, H)
      @test index(A, AA) == 2
      @test is_normal_subgroup(AA, A)
      @test is_normalized_by(AA, A)
      @test phi(AA[1]) == AA[1]
      @test phi(AA[2]) == AA[2]
      @test order(quo(A, AA)[1]) == 2
      @test is_invariant(f, H)

      S = sylow_subgroup(G, 3)[1]
      x = gen(S, 1)
      f = A(hom(G, G, y -> y^x))
      fHa = restrict_automorphism(f, H)
      fHh = restrict_homomorphism(f, H)
      @test parent(fHa) == automorphism_group(H)
      @testset for g in gens(H)
         @test fHa(g) == H(f(g))
         @test fHh(g) == H(f(g))
      end

      V, _ = pcore(G, 2)
      S, g = quo(G, V)
      @test induced_automorphism(g, f) == automorphism_group(S)(inner_automorphism(g(x)))
   end
end

@testset "Other automorphisms groups" begin
   C = cyclic_group(3)
   G = direct_product(C,C)
   A = automorphism_group(G)

   @test is_isomorphic(A,GL(2,3))
   @test order(inner_automorphism_group(A)[1])==1

   # Create an Oscar group from a group of automorphisms in GAP.
   G = alternating_group(6)
   A = automorphism_group(G)
   B = Oscar._oscar_group(GapObj(A))
   @test B == A
   @test B !== A
   @test B.X === A.X

   F = free_group(2)
   x, y = gens(F)
   Q, = quo(F, [x^-3*y^-3*x^-1*y*x^-2,
                x*y^-1*x^-1*y*x^-1*y^-1*x^3*y^-1,
                x*y^-1*x^-1*y^2*x^-2*y^-1*x^2])
   A = automorphism_group(Q)
   q = gen(Q, 1)
   @test all(a -> a(preimage(a, q)) == q, collect(A))
end

@testset "Composition of mappings" begin
   g = symmetric_group(4)
   q, epi = quo(g, pcore(g, 2)[1])
   iso = @inferred isomorphism(PermGroup, q)
   comp = compose(epi, iso)
   @test domain(comp) == domain(epi)
   @test codomain(comp) == codomain(iso)
end
