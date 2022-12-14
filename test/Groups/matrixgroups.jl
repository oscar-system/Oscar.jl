@testset "Oscar-GAP relationship for finite fields" begin
   F = GF(29, 1)
   z = F(2)
   G = GL(3,F)
   @test G.X isa GAP.GapObj
   @test isdefined(G,:X)
   @test isdefined(G, :ring_iso)
   @test G.ring_iso(z) isa GAP.FFE
   Z = G.ring_iso(z)
   @test Z in codomain(G.ring_iso)
   @test preimage(G.ring_iso, Z)==z
   @test domain(G.ring_iso) == F
   @test GAP.Globals.IsField(codomain(G.ring_iso))
   @test GAP.Globals.Size(codomain(G.ring_iso))==29
   @test GAP.Globals.IsZero(14*Z+1)
   @test iszero(preimage(G.ring_iso, GAP.Globals.Zero(codomain(G.ring_iso))))
   @test GAP.Globals.IsOne(G.ring_iso(one(F)))
   @test isone(preimage(G.ring_iso, GAP.Globals.One(codomain(G.ring_iso))))
   
   xo = matrix(F,3,3,[1,z,0,0,1,2*z+1,0,0,z+2])
#   xg = Vector{GAP.GapObj}(undef, 3)
#   for i in 1:3
#      xg[i] = GAP.GapObj([preimage(G.ring_iso, xo[i,j]) for j in 1:3])
#   end
#   xg=GAP.Obj(xg)

   xg = GAP.GapObj([[G.ring_iso(xo[i,j]) for j in 1:3] for i in 1:3]; recursive=true)
   @test map_entries(G.ring_iso, xo) == xg
   @test Oscar.preimage_matrix(G.ring_iso, xg) == xo
   @test Oscar.preimage_matrix(G.ring_iso, GAP.Globals.One(GAP.Globals.GL(3, codomain(G.ring_iso)))) == matrix(one(G))
   @test GAP.Globals.Order(map_entries(G.ring_iso, diagonal_matrix([z,z,one(F)]))) == 28

   T,t = PolynomialRing(GF(3) ,"t")
   F,z = FiniteField(t^2+1,"z")
   G = GL(3,F)
   @test G.X isa GAP.GapObj
   @test isdefined(G,:X)
   @test isdefined(G, :ring_iso)
   @test G.ring_iso(z) isa GAP.FFE
   Z = G.ring_iso(z)
   @testset for a in F
      for b in F
         @test G.ring_iso(a*b)==G.ring_iso(a)*G.ring_iso(b)
         @test G.ring_iso(a-b)==G.ring_iso(a)-G.ring_iso(b)
      end
   end
   @test Z in codomain(G.ring_iso)
   @test preimage(G.ring_iso, Z)==z
   @test preimage(G.ring_iso, G.ring_iso(F(2)))==F(2)
   @test domain(G.ring_iso) == F
   @test GAP.Globals.IsField(codomain(G.ring_iso))
   @test GAP.Globals.Size(codomain(G.ring_iso))==9
   @test iszero(preimage(G.ring_iso, GAP.Globals.Zero(codomain(G.ring_iso))))
   @test GAP.Globals.IsZero(Z^2+1)
   @test GAP.Globals.IsOne(G.ring_iso(one(F)))
   @test isone(preimage(G.ring_iso, GAP.Globals.One(codomain(G.ring_iso))))
   
   xo = matrix(F,3,3,[1,z,0,0,1,2*z+1,0,0,z+2])
   xg = Vector{GAP.GapObj}(undef, 3)
   for i in 1:3
      xg[i] = GAP.Obj([G.ring_iso(xo[i,j]) for j in 1:3])
   end
   xg=GAP.Obj(xg)
   @test map_entries(G.ring_iso, xo) == xg
   @test Oscar.preimage_matrix(G.ring_iso, xg) == xo
   @test Oscar.preimage_matrix(G.ring_iso, GAP.Globals.One(GAP.Globals.GL(3, codomain(G.ring_iso)))) == matrix(one(G))
   @test GAP.Globals.Order(map_entries(G.ring_iso, diagonal_matrix([z,z,one(F)])))==4
end

@testset "Oscar-GAP relationship for cyclotomic fields" begin
   fields = Any[CyclotomicField(n) for n in [1, 3, 4, 5, 8, 15, 45]]
   push!(fields, (QQ, QQ(1)))
   F, z = abelian_closure(QQ)
   push!(fields, (F, z(5)))

   @testset for (F, z) in fields
      f = Oscar.iso_oscar_gap(F)
      g = elm -> map_entries(f, elm)
      G = MatrixGroup(3, F)
      mats = [matrix(F, [0 z 0; 0 0 1; 1 0 0]),
              matrix(F, [0 1 0; 1 0 0; 0 0 1])]
      G.gens = [MatrixGroupElem(G, m) for m in mats]
      for a in map(matrix, gens(G))
         for b in map(matrix, gens(G))
            @test g(a * b) == g(a) * g(b)
            @test g(a - b) == g(a) - g(b)
         end
      end
      @test G.ring_iso(z) isa GAP.Obj
      @test G.X isa GAP.GapObj
      @test isdefined(G, :X)
      @test isdefined(G, :ring_iso)
      Z = G.ring_iso(z)
      @test Z in codomain(G.ring_iso)
      @test preimage(G.ring_iso, Z) == z
      @test domain(G.ring_iso) == F
      @test GAP.Globals.IsField(codomain(G.ring_iso))
      @test iszero(preimage(G.ring_iso, GAP.Globals.Zero(codomain(G.ring_iso))))
      @test GAP.Globals.IsOne(G.ring_iso(one(F)))
      @test isone(preimage(G.ring_iso, GAP.Globals.One(codomain(G.ring_iso))))

      xo = matrix(F, 3, 3, [0, 1, 0, 0, 1, z, 0, 0, z])
      xg = GAP.GapObj([[G.ring_iso(xo[i, j]) for j in 1:3] for i in 1:3]; recursive = true)
      @test map_entries(G.ring_iso, xo) == xg
      @test Oscar.preimage_matrix(G.ring_iso, xg) == xo
      @test Oscar.preimage_matrix(G.ring_iso, GAP.Globals.IdentityMat(3)) == matrix(one(G))
      if F isa AnticNumberField
         flag, n = Hecke.is_cyclotomic_type(F)
         @test GAP.Globals.Order(map_entries(G.ring_iso, diagonal_matrix([z, z, one(F)]))) == n
      end
   end
end

@testset "faithful reduction from char. zero to finite fields" begin

   M = matrix(QQ, [ 2 0; 0 2 ])
   @test_throws ErrorException Oscar.isomorphic_group_over_finite_field(matrix_group([M]))

   K, a = CyclotomicField(5, "a")
   L, b = CyclotomicField(3, "b")

   inputs = [
     #[ matrix(ZZ, [ 0 1 0; -1 0 0; 0 0 -1 ]) ],
     [ matrix(QQ, [ 0 1 0; -1 0 0; 0 0 -1 ]) ],
     [ matrix(K, [ a 0; 0 a ]) ],
     [ matrix(L, 2, 2, [ b, 0, -b - 1, 1 ]), matrix(L, 2, 2, [ 1, b + 1, 0, b ]) ]
   ]

   @testset "... over ring $(base_ring(mats[1]))" for mats in inputs
     G0 = matrix_group(mats)
     G, g = Oscar.isomorphic_group_over_finite_field(G0)

     @test !has_order(G0)
     order(G0)
     @test has_order(G0)

     for i in 1:10
       x, y = rand(G), rand(G)
       @test (g\x) * (g\y) == g\(x * y)
       @test g(g\x) == x
     end

     H = GAP.Globals.Group(GAP.Obj(gens(G0); recursive=true))
     f = GAP.Globals.GroupHomomorphismByImages(G.X, H)
     @test GAP.Globals.IsBijective(f)
     @test order(G) == GAP.Globals.Order(H)
   end

   G = MatrixGroup(2, QQ, dense_matrix_type(QQ)[])
   @test order(Oscar.isomorphic_group_over_finite_field(G)[1]) == 1
end

@testset "Type operations" begin
   G = GL(5,5)
   x = rand(G)
   @test ring_elem_type(typeof(G))==typeof(one(base_ring(G)))
   @test mat_elem_type(typeof(G))==typeof(matrix(x))
   @test elem_type(typeof(G))==typeof(x)
   @test Oscar._gap_filter(typeof(G))(G.X)
end

#FIXME : this may change in future. It can be easily skipped.
@testset "Fields assignment" begin
   T,t=PolynomialRing(GF(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   @test G isa MatrixGroup
   @test F==base_ring(G)
   @test 2==degree(G)
   @test !isdefined(G,:X)
   @test !isdefined(G,:gens)
   @test !isdefined(G,:ring_iso)

   @test order(G)==5760
   @test order(G) isa fmpz
   @test order(Int64, G) isa Int64
   @test order(Int64, SL(2,F))==720
   @test isdefined(G,:X)
   @test !isdefined(G,:gens)

   @test ngens(G)==2
   @test isdefined(G,:gens)
   @test typeof(gens(G)) == Vector{elem_type(G)}

   x = matrix(F,2,2,[1,0,0,1])
   x = G(x)
   @test !isdefined(x,:X)
   @test x.X isa GAP.GapObj
   x = matrix(G[1])
   x = G(x)
   @test !isdefined(x,:X)
   x = matrix(F,2,2,[1,z,z,1])
   x = G(x)
   @test !isdefined(x,:X)
   @test order(x)==8
   @test isdefined(x,:X)

   x = G([2,1,2,0])
   y = G([z+1,0,0,z+2])
   @test parent(x)==G
   H,f = sub(G,[x,y])
   @test isdefined(H,:gens)
   @test gens(H)==[x,y]
   @test typeof(gens(H)) == Vector{elem_type(H)}
   @test H==SL(2,F)
   @test parent(x)==G
   @test parent(H[1])==H
   @test parent(f(H[1]))==G

   K1 = matrix_group(x,y,x*y)
   @test K1.X isa GAP.GapObj
   @test K1.X==H.X

   K = matrix_group(x,x^2,y)
   @test isdefined(K, :gens)
   @test !isdefined(K,:X)
   @test K.gens==[x,x^2,y]
   @test typeof(gens(K)) == Vector{elem_type(K)}
   @test parent(x)==G
   @test x==K[1]                           #TODO changes in future if we decide to keep track of the parent
   @test parent(K[1])==K
   @test H==K
   @test H.gens != K.gens
   @test K==matrix_group([x,x^2,y])
   @test K==matrix_group(matrix(x), matrix(x^2), matrix(y))
   @test K==matrix_group([matrix(x), matrix(x^2), matrix(y)])

   G = MatrixGroup(nrows(x), F)
   G.gens = typeof(x)[]   # empty list of generators
   G.X
   @test one(G) == one(x)

   G = GL(3,F)
   x = G([1,z,0,0,z,0,0,0,z+1])
   @test order(x)==8

   G = MatrixGroup(4,F)
   @test_throws ErrorException G.X
   setfield!(G,:descr,:GX)
   @test isdefined(G,:descr)
   @test_throws ErrorException G.X
end


@testset "Constructors" begin
   @testset for n in 4:5
      @testset for F in [GF(2, 2), GF(3, 1)]
         q = Int(order(F))
         G = GL(n,F)
         S = SL(n,F)
         @test G==GL(n,q)
         @test G==general_linear_group(n,F)
         @test G==general_linear_group(n,q)
         @test S==SL(n,q)
         @test S==special_linear_group(n,F)
         @test S==special_linear_group(n,q)
         @test order(S)==prod(BigInt[q^n-q^i for i in 0:(n-1)])÷(q-1)
         @test index(G,S)==q-1
      end
   end

   @testset for n in 1:3
      @testset for q in [2,3,4]
         @test unitary_group(n,q)==GU(n,q)
         @test special_unitary_group(n,q)==SU(n,q)
         @test index(GU(n,q),SU(n,q))==q+1
      end
   end

   @testset for n in [4,6]
      @testset for F in [GF(2, 2), GF(3, 1)]
         q = Int(order(F))
         G = Sp(n,F)
         @test G==Sp(n,q)
         @test G==symplectic_group(n,F)
         @test G==symplectic_group(n,q)
      end
   end

   @testset for F in [GF(3, 1), GF(2, 2), GF(5, 1)]
      q = Int(order(F))
      @testset for n in [4,6]
         @testset for e in [+1,-1]
            G = GO(e,n,F)
            S = SO(e,n,F)
            O = omega_group(e,n,F)
            @test G==GO(e,n,q)
            @test G==orthogonal_group(e,n,F)
            @test G==orthogonal_group(e,n,q)
            @test S==SO(e,n,q)
            @test S==special_orthogonal_group(e,n,F)
            @test S==special_orthogonal_group(e,n,q)
            @test O==omega_group(e,n,q)
            @test index(S,O)==2
            @test index(G,S) == gcd(2, q-1)
         end
      end
      @testset for n in [3,5]
         G = GO(n,F)
         S = SO(n,F)
         O = omega_group(n,F)
         @test G==GO(n,q)
         @test G==orthogonal_group(n,F)
         @test G==orthogonal_group(n,q)
         @test S==SO(n,q)
         @test S==special_orthogonal_group(n,F)
         @test S==special_orthogonal_group(n,q)
         @test O==omega_group(n,q)
         @test index(G,S) == gcd(2, q-1)
         @test index(S,O) == gcd(2, q-1)
      end
   end

   @test order(omega_group(+1,4,3))==288
   @test order(omega_group(-1,4,3))==360
   @test order(omega_group(3,3))==12
   @test order(omega_group(3,4)) == 60

   @test_throws ArgumentError GO(+2,2,5)
   @test_throws ArgumentError SO(+2,2,5)
   @test_throws ArgumentError omega_group(-2,4,3)

   @test omega_group(1,5)==SO(1,5)
   @test index(GO(1,7),omega_group(1,7))==2
   @test order(omega_group(1,5))==1
   G = omega_group(1,4,2)
   @testset for x in gens(G)
       @test iseven(rank(matrix(x)-1))
   end
end


@testset "Assignments and generators" begin

   G = GL(3,5)
   S = SL(3,5)
   x1 = G([-1,0,1,-1,0,0,0,-1,0])
   x2 = G([2,0,0,0,3,0,0,0,1])
   @test x1==G([4,0,1,4,0,0,0,4,0])
   @test x1==G([4 0 1; 4 0 0; 0 4 0])
   @test matrix_group(x1,x2)==MatrixGroup(3,base_ring(x1),[x1,x2])
   @test matrix_group(x1,x2)==matrix_group([x1,x2])
   H = matrix_group([x1,x2])
   @test isdefined(H,:gens)
   @test H[1]==x1
   @test H[2]==x2
   @test parent(H[1])==H
   @test parent(x1)==G
   @test H==S
   H1 = matrix_group([x1,x2])
   @test H==H1
   @test !isdefined(H1,:X)
   H1 = matrix_group([matrix(x1),matrix(x2)])
   @test H==H1
   @test parent(H1[1])==H1
   @test !isdefined(H1,:X)
   H1 = matrix_group(x1,x2)
   @test H==H1
   @test parent(H1[1])==H1
   @test !isdefined(H1,:X)
   H1 = matrix_group(matrix(x1),matrix(x2))
   @test H==H1
   @test parent(H1[1])==H1
   @test !isdefined(H1,:X)
   x3 = matrix(base_ring(G),3,3,[0,0,0,0,1,0,0,0,1])
   @test_throws AssertionError matrix_group(matrix(x1),x3)
   @test parent(x1)==G

   G4 = GL(4,5)
   x3 = G4([1,0,2,0,0,1,0,2,0,0,1,0,0,0,0,1])
   @test_throws AssertionError matrix_group([x1,x3])

   G4 = GL(3,7)
   x3 = G4([2,0,0,0,3,0,0,0,1])
   @test_throws AssertionError matrix_group([x1,x3])
end

@testset "Iterator" begin
   G = SL(2,3)
   N = 0
   for x in G
      N+=order(x)
   end
   @test N==99

   @test Set(collect(G))==Set([x for x in G])
end

@testset "Membership" begin
   T,t=PolynomialRing(GF(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   S = SL(2,F)
   O = GO(1,2,F)

   x = matrix(F,2,2,[1,z,0,z])
   @test x in G
   @test !(x in S)
   @test !(matrix(F,2,2,[0,0,0,1]) in G)
   @test_throws ArgumentError S(x)
   @test G(x) isa MatrixGroupElem
   @test S(x; check=false)==G(x)
   @test S(G(x); check=false)==G(x)
   x = G(x)
   y = MatrixGroupElem(G,x.X)
   @test_throws ArgumentError S(y)
   @test G(y) isa MatrixGroupElem
   @test G(y*y)==G(y)*G(y)
   @test x==G([1,z,0,z])
   @test x==G([1 z; 0 z])
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

   x = matrix(O[1]*O[2])
   @test x in G
   @test x in O
   @test parent(O(x))==O
   x = G(x)
   @test parent(O(x))==O
   @test_throws ArgumentError O([z,0,0,1])

   K = matrix_group(S[1],S[2],S[1]*S[2])
   x = S(matrix(F,2,2,[2,z,0,2]))
   @test x in K
   @test isdefined(K,:X)
   @test isdefined(x,:X)

end

@testset "Methods on elements" begin
   T,t=PolynomialRing(GF(3),"t")
   F,z=FiniteField(t^2+1,"z")

   G = GL(2,F)
   x = G([1,z,0,1])
   y = G([z+2,0,0,1])
   @test y[1,1]==z+2
   @test x[1,2]==z
   @test x*y==G([z+2,z,0,1])
   @test x^-1==G([1,2*z,0,1])
   @test y^x==G([z+2,z+2,0,1])
   @test comm(y,x)==G([1,1,0,1])
   @test isone(G([1,0,0,1]))
   @test !isone(x)
   @test det(x)==1
   @test det(y)==z+2
   @test tr(x)==2
   @test tr(y)==z
   @test order(x)==3
   @test order(y)==8
   @test base_ring(x)==F
   @test nrows(y)==2
   @test x*matrix(y) isa fq_nmod_mat
   @test matrix(x*y)==matrix(x)*y
   @test G(x*matrix(y))==x*y
   @test matrix(x)==x.elm

   xg = GAP.Globals.Random(G.X)
   yg = GAP.Globals.Random(G.X)
   pg = MatrixGroupElem(G, xg*yg)
   @test pg == MatrixGroupElem(G, Oscar.preimage_matrix(G.ring_iso, xg))*MatrixGroupElem(G, Oscar.preimage_matrix(G.ring_iso, yg))

   O = GO(-1,2,F)
   S = SL(2,F)
   xo = O([2,0,2,1])
   xs = S([1,z,0,1])
   @test parent(xs*xo)==G
   xo = xo^2
   @test parent(xs*xo)==G
   @test parent(xs*S(xo))==S
end

@testset "Subgroups" begin
   T,t=PolynomialRing(GF(3),"t")
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
   @test is_normal(O,H)
   @test index(O,H)==2
#   @test index(GO(0,3,3), omega_group(0,3,3))==4
   @test index(GO(1,2,8), omega_group(1,2,8))==2
   @test index(GO(1,2,9), omega_group(1,2,9))==4
end

@testset "Cosets and conjugacy classes" begin
   T,t=PolynomialRing(GF(3),"t")
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
   @test x^G[2] in collect(cc)
   @test representative(cc)==x
   @test parent(representative(cc))==G
   @test length(cc)==index(G,C)

   cc = conjugacy_class(G,H)
   @test H^G[2] in collect(cc)
   @test representative(cc)==H
   @test length(cc)==index(G,normalizer(G,H)[1])
   @test rand(cc) in collect(cc)

   x = G([1,z,0,1])
   y = G([1,0,0,z+1])
   H = matrix_group([x])
   @test gens(H)==[x]
   K = H^y
   @test gens(K)==[x^y]
   @test !isdefined(K,:X)


   G = GL(2,3)
   @test length(conjugacy_classes(G))==8
   @test length(@inferred conjugacy_classes_subgroups(G))==16
   @test length(@inferred conjugacy_classes_maximal_subgroups(G))==3
end

@testset "Jordan structure" begin
   F = GF(3, 1)
   R,t = PolynomialRing(F,"t")
   G = GL(9,F)

   L_big = [
        [(t-1,3), (t^2+1,1), (t^2+1,2)],
        [(t^9 + t^7 + 2*t^6 + t^5 + 2*t^4 + t^3 + 2*t^2 + 2*t + 1, 1)],
        [(t-2,9)]
       ]
   @testset for L in L_big
      x = cat([generalized_jordan_block(a...) for a in L]..., dims=(1,2))
   # TODO: the change_base_ring is necessary, otherwise the equality between polynomials does not work
      @test MSet([(change_base_ring(F,f[1]),f[2]) for f in pol_elementary_divisors(x) ])==MSet([(change_base_ring(F,f[1]),f[2]) for f in L])
      @test MSet([(change_base_ring(F,f[1]),f[2]) for f in pol_elementary_divisors(G(x)) ])==MSet([(change_base_ring(F,f[1]),f[2]) for f in L])
      s, u = multiplicative_jordan_decomposition(G(x))
      @test parent(s)==G
      @test parent(u)==G
      @test is_coprime(order(s),3)
      @test isone(u) || is_power(order(u))[2]==3
      @test is_semisimple(s)
      @test is_unipotent(u)
      @test s*u==G(x)
      @test s*u==u*s

      z = matrix(rand(G))
      x = z^-1*x*z
      a,b = generalized_jordan_form(x)
      @test b^-1*a*b==x
      z = matrix(rand(G))
      @test generalized_jordan_form(z^-1*x*z)[1]==a
      @test generalized_jordan_form(a)[1]==a
   end

   x = one(G)
   @test is_semisimple(x) && is_unipotent(x)

   F,z = FiniteField(5,3,"z")
   G = GL(6,F)
   R,t = PolynomialRing(F,"t")
   f = t^3+t*z+1
   x = generalized_jordan_block(f,2)
   @test generalized_jordan_block(f,2)==hvcat((2,2),companion_matrix(f),identity_matrix(F,3),zero_matrix(F,3,3),companion_matrix(f))
   @testset for i in [2,4,42,62]
      y = Oscar._elem_given_det(G(x),z^i)
      @test x*y==y*x
      @test det(y)==z^i
   end
   @test_throws ErrorException Oscar._elem_given_det(G(x),z)

   @testset "Low-level methods in linear_centralizer.jl" begin
      @test Oscar._SL_order(3,fmpz(8))== fmpz(div(prod([8^3-8^i for i in 0:2]),7))
      @test Oscar._SL_order(4, GF(3, 1))== fmpz(div(prod([3^4-3^i for i in 0:3]),2))
      L = Oscar._gens_for_GL(1,GF(7, 1))
      @test length(L)==1
      @test L[1]^2 !=1 && L[1]^3 !=1
      L = Oscar._gens_for_GL(4,GF(2, 2))
      @test length(L)==2
      @test matrix_group(L...)==GL(4,GF(2, 2))
      L = Oscar._gens_for_SL(5,GF(3, 1))
      @test matrix_group(L...)==SL(5,GF(3, 1))
      L = Oscar._gens_for_GL(5,GF(2, 1))
      @test length(L)==2
      @test matrix_group(L...)==GL(5,GF(2, 1))
      _,t = PolynomialRing(GF(3, 1),"t")
      f = t^2+t-1
      L = Oscar._gens_for_GL_matrix(f,2,GF(3, 1); D=2)
      @test length(L)==2
      @test nrows(L[1])==8
      @test L[1]^8==1
      @test L[2]^3==1
      @test order(matrix_group(L...))==order(GL(2,9))
      L = Oscar._gens_for_SL_matrix(f,2,GF(3, 1); D=2)
      @test length(L)==3
      @test nrows(L[1])==8
      @test L[1]^8==1
      @test L[2]^3==1
      @test order(matrix_group(L...))==div(order(GL(2,9)),2)
      x = cat([generalized_jordan_block(f,n) for n in [1,1,1,2,2,3]]..., dims=(1,2))
      L,c = Oscar._centr_block_unipotent(f,GF(3, 1),[1,1,1,2,2,3])
      @testset for l in L
         @test l*x==x*l
      end
      @test c==order(GL(3,9))*order(GL(2,9))*8*BigInt(9)^32
   end
end

@testset "isometry group" begin
   q = quadratic_space(QQ,QQ[2 1; 1 2])
   L = lattice(q, QQ[1 0; 0 1])
   G = isometry_group(L)
   @test order(G) == 12
   @test isometry_group(L) == orthogonal_group(L)
   # L = lattice(q, QQ[0 0; 0 0], isbasis=false)
   # @test order(isometry_group(L)) == 1

   Qx, x = PolynomialRing(FlintQQ, "x", cached = false)
   f = x^2-2;
   K, a = number_field(f)
   D = matrix(K, 3, 3, [2, 0, 0, 0, 1, 0, 0, 0, 7436]);
   gens = [[13, 0, 0], [156*a+143, 0, 0], [3//2*a+5, 1, 0], [3//2*a+5, 1, 0], [21//2*a, 0, 1//26], [21//2*a, 0, 1//26]]
   L = quadratic_lattice(K, gens, gram = D)
   G = orthogonal_group(L)
   g = -identity_matrix(K, 3)
   @test g in G
end

@testset "deepcopy" begin
   g = general_linear_group(2, 4)

   m = MatrixGroupElem(g, gen(g, 1).X);  # do not call `show`!
   @test isdefined(m, :X)
   @test ! isdefined(m, :elm)
   c = deepcopy(m);
   @test isdefined(c, :X)
   @test ! isdefined(c, :elm)
   @test c.X == m.X

   m = MatrixGroupElem(g, matrix(gen(g, 1)), gen(g, 1).X)
   @test isdefined(m, :X)
   @test isdefined(m, :elm)
   c = deepcopy(m);
   @test isdefined(c, :X)
   @test isdefined(c, :elm)
   @test c.X == m.X
   @test matrix(c) == matrix(m)

   m = MatrixGroupElem(g, matrix(gen(g, 1)))
   @test ! isdefined(m, :X)
   @test isdefined(m, :elm)
   c = deepcopy(m);
   @test ! isdefined(c, :X)
   @test isdefined(c, :elm)
   @test matrix(c) == matrix(m)

   @test deepcopy([one(g)]) == [one(g)]
end

@testset "matrix action on vectors" begin
   for R in [ZZ, QQ, GF(2,2)]
     T = elem_type(R)
     mat = matrix(R, [1 0; 0 1])
     G = matrix_group([mat])
     h = gens(G)[1]
     v = [R(x) for x in [1, 1]]
     @test v * h isa Vector{T}
     @test v * h == v * mat
     @test h * v isa Vector{T}
     @test h * v == mat * v
   end
end
