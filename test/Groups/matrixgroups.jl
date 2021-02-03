@testset "Oscar-GAP relationship" begin
   F = GF(29,1)[1]
   z = F(2)
   G = GL(3,F)
   @test G.X isa GapObj
   @test isdefined(G,:X)
   @test isdefined(G, :ring_iso)
   @test isdefined(G, :mat_iso)
   @test G.ring_iso(z) isa FFE
   Z = G.ring_iso(z)
   @test GAP.Globals.IN(Z,G.ring_iso.codomain)
   @test G.ring_iso(Z)==z
   @test G.ring_iso.domain==F
   @test GAP.Globals.IsField(G.ring_iso.codomain)
   @test GAP.Globals.Size(G.ring_iso.codomain)==29
   @test GAP.Globals.IsZero(14*Z+1)
   @test iszero(G.ring_iso(GAP.Globals.Zero(G.ring_iso.codomain)))
   @test GAP.Globals.IsOne(G.ring_iso(one(F)))
   @test isone(G.ring_iso(GAP.Globals.One(G.ring_iso.codomain)))
   
   xo = matrix(F,3,3,[1,z,0,0,1,2*z+1,0,0,z+2])
#   xg = Vector{GapObj}(undef, 3)
#   for i in 1:3
#      xg[i] = GapObj([G.ring_iso(xo[i,j]) for j in 1:3])
#   end
#   xg=GAP.julia_to_gap(xg)
   riso=Oscar.gen_ring_iso(F)
   @test riso isa Oscar.GenRingIso
   xg = Oscar.mat_oscar_gap(xo,riso)
   @test xg isa GapObj
   @test Oscar.mat_gap_oscar(xg,riso)==xo
   @test Oscar.mat_oscar_gap(Oscar.mat_gap_oscar(xg,riso))==xg

   xg = GapObj([[G.ring_iso(xo[i,j]) for j in 1:3] for i in 1:3]; recursive=true)
   @test G.mat_iso(xo) isa GapObj
   @test G.mat_iso(xo)==xg
   @test G.mat_iso(xg)==xo
   @test G.mat_iso(GAP.Globals.One(GAP.Globals.GL(3,G.ring_iso.codomain)))==one(G).elm
   @test GAP.Globals.Order(G.mat_iso(diagonal_matrix([z,z,one(F)])))==28

   T,t = PolynomialRing(GF(3),"t")
   F,z = FiniteField(t^2+1,"z")
   G = GL(3,F)
   @test G.X isa GapObj
   @test isdefined(G,:X)
   @test isdefined(G, :ring_iso)
   @test isdefined(G, :mat_iso)
   @test G.ring_iso(z) isa FFE
   Z = G.ring_iso(z)
   @test GAP.Globals.IN(Z,G.ring_iso.codomain)
   @test G.ring_iso(Z)==z
   @test G.ring_iso(G.ring_iso(F(2)))==F(2)
   @test G.ring_iso.domain==F
   @test GAP.Globals.IsField(G.ring_iso.codomain)
   @test GAP.Globals.Size(G.ring_iso.codomain)==9
   @test iszero(G.ring_iso(GAP.Globals.Zero(G.ring_iso.codomain)))
   @test GAP.Globals.IsZero(Z^2+1)
   @test GAP.Globals.IsOne(G.ring_iso(one(F)))
   @test isone(G.ring_iso(GAP.Globals.One(G.ring_iso.codomain)))
   
   xo = matrix(F,3,3,[1,z,0,0,1,2*z+1,0,0,z+2])
   xg = Vector{GapObj}(undef, 3)
   for i in 1:3
      xg[i] = GAP.julia_to_gap([G.ring_iso(xo[i,j]) for j in 1:3])
   end
   xg=GAP.julia_to_gap(xg)
   @test G.mat_iso(xo)==xg
   @test G.mat_iso(xg)==xo
   @test G.mat_iso(GAP.Globals.One(GAP.Globals.GL(3,G.ring_iso.codomain)))==one(G).elm
   @test GAP.Globals.Order(G.mat_iso(diagonal_matrix([z,z,one(F)])))==4
end

@testset "Type operations" begin
   G = GL(5,5)
   x = rand(G)
   @test ring_elem_type(typeof(G))==typeof(one(base_ring(G)))
   @test mat_elem_type(typeof(G))==typeof(x.elm)
   @test elem_type(typeof(G))==typeof(x)
   @test Oscar._gap_filter(typeof(G))(G.X)
end

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
   @test !isdefined(G,:mat_iso)

   @test order(G)==5760
   @test isdefined(G,:X)
   @test !isdefined(G,:gens)

   @test ngens(G)==2
   @test isdefined(G,:gens)

   x = matrix(F,2,2,[1,0,0,1])
   x = G(x)
   @test !isdefined(x,:X)
   @test x.X isa GapObj
   x = G[1].elm
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
   @test H==SL(2,F)
   @test parent(x)==G
   @test parent(H[1])==H
   @test parent(f(H[1]))==G

   K1 = matrix_group(x,y,x*y)
   @test K1.X isa GapObj
   @test K1.X==H.X

   K = matrix_group(x,x^2,y)
   @test isdefined(K, :gens)
   @test !isdefined(K,:X)
   @test K.gens==[x,x^2,y]
   @test parent(x)==G
   @test x==K[1]                           #TODO changes in future if we decide to keep track of the parent
   @test parent(K[1])==K
   @test H==K
   @test H.gens != K.gens
   @test K==matrix_group([x,x^2,y])
   @test K==matrix_group(x.elm, (x^2).elm, y.elm)
   @test K==matrix_group([x.elm, (x^2).elm, y.elm])


   G = GL(3,F)
   x = G([1,z,0,0,z,0,0,0,z+1])
   @test order(x)==8
   @test isdefined(G,:mat_iso)

   G = GL(4,2)
   @test G.mat_iso isa Oscar.GenMatIso
   
   G = MatrixGroup(4,F)
   @test_throws ErrorException G.X
   setfield!(G,:descr,:GX)
   @test isdefined(G,:descr)
   @test_throws ErrorException G.X
end


@testset "Constructors" begin
   @testset for n in 4:5
      @testset for F in [GF(2,2)[1], GF(3,1)[1]]
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
      @testset for F in [GF(2,2)[1], GF(3,1)[1]]
         q = Int(order(F))
         G = Sp(n,F)
         @test G==Sp(n,q)
         @test G==symplectic_group(n,F)
         @test G==symplectic_group(n,q)
      end
   end

   @testset for F in [GF(3,1)[1], GF(2,2)[1], GF(5,1)[1]]
      q = Int(order(F))
      @testset for n in [4,6]
         @testset for e in [+1,-1]
            G = GO(e,n,F)
            S = SO(e,n,F)
            @test G==GO(e,n,q)
            @test G==orthogonal_group(e,n,F)
            @test G==orthogonal_group(e,n,q)
            @test S==SO(e,n,q)
            @test S==special_orthogonal_group(e,n,F)
            @test S==special_orthogonal_group(e,n,q)
            if isodd(q)
               @test index(G,S)==2
               @test order(S)==2*order(omega_group(e,n,q))
            else
               @test order(G)==2*order(omega_group(e,n,q))
            end
         end
      end
      @testset for n in [3,5]
         if isodd(q)
            G = GO(n,F)
            S = SO(n,F)
            @test G==GO(n,q)
            @test G==orthogonal_group(n,F)
            @test G==orthogonal_group(n,q)
            @test S==SO(n,q)
            @test S==special_orthogonal_group(n,F)
            @test S==special_orthogonal_group(n,q)
            @test index(G,S)==2
            @test order(S)==2*order(omega_group(n,q))
         end
      end
   end

   @test order(omega_group(+1,4,3))==288
   @test order(omega_group(-1,4,3))==360
   @test order(omega_group(3,3))==12

   G = GL(4,3)
   
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
   H1 = matrix_group([x1.elm,x2.elm])
   @test H==H1
   @test parent(H1[1])==H1
   @test !isdefined(H1,:X)
   H1 = matrix_group(x1,x2)
   @test H==H1
   @test parent(H1[1])==H1
   @test !isdefined(H1,:X)
   H1 = matrix_group(x1.elm,x2.elm)
   @test H==H1
   @test parent(H1[1])==H1
   @test !isdefined(H1,:X)
   x3 = matrix(base_ring(G),3,3,[0,0,0,0,1,0,0,0,1])
   @test_throws AssertionError matrix_group(x1.elm,x3)
   @test parent(x1)==G

   G4 = GL(4,5)
   x3 = G4([1,0,2,0,0,1,0,2,0,0,1,0,0,0,0,1])
   @test_throws AssertionError matrix_group([x1,x3])

   G4 = GL(3,7)
   x3 = G4([2,0,0,0,3,0,0,0,1])
   @test_throws AssertionError matrix_group([x1,x3])
end


@testset "Classical groups as PermGroup" begin
   @test GL(PermGroup,2,3) isa PermGroup
   @test order(G) isa fmpz
   @test order(Int64, G) isa Int64
   @test order(GL(PermGroup,2,3)) == order(GL(2,3))
   @test isisomorphic(GL(PermGroup,2,3),GL(2,3))[1]
   F2 = GF(2,1)[1]
   F3 = GF(3,1)[1]
   @test SL(PermGroup,2,F2) isa PermGroup

   @test order(SL(PermGroup,2,2)) == order(SL(2,2))
   @test order(SL(PermGroup,2,F2)) == order(SL(2,2))
   @test order(GU(PermGroup,2,2)) == order(GU(2,2))
   @test order(SU(PermGroup,2,2)) == order(SU(2,2))
   @test order(GO(PermGroup,+1,2,2)) == order(GO(+1,2,2))
   @test order(SO(PermGroup,+1,2,3)) == order(SO(+1,2,3))
   @test order(GO(PermGroup,+1,2,F2)) == order(GO(+1,2,2))
   @test order(SO(PermGroup,+1,2,F3)) == order(SO(+1,2,3))
   @test order(Sp(PermGroup,2,2)) == order(Sp(2,2))
   @test SL(MatrixGroup,2,3) == SL(2,3)
   @test GL(MatrixGroup,2,3) == GL(2,3)
   @test Sp(MatrixGroup,2,3) == Sp(2,3)
   @test GO(MatrixGroup,3,3) == GO(3,3)
   @test GO(MatrixGroup,-1,2,3) == GO(-1,2,3)
   @test SO(MatrixGroup,-1,2,3) == SO(-1,2,3)
   @test GU(MatrixGroup,2,3) == GU(2,3)
   @test SU(MatrixGroup,2,3) == SU(2,3)
end


@testset "Iterator" begin
   G = SL(2,3)
   N = 0
   for x in G
      N+=order(x)
   end
   @test N==99

   @test Set(elements(G))==Set([x for x in G])
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
   @test !(matrix(F,2,2,[0,0,0,1]) in G)
   @test_throws ArgumentError S(x)
   @test G(x) isa MatrixGroupElem
   @test S(x; check=false)==G(x)
   x = G(x)
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

   x = (O[1]*O[2]).elm
   @test x in G
   @test x in O
   @test parent(O(x))==O
   x = G(x)
   @test parent(O(x))==O
   @test_throws ArgumentError O([z,0,0,1])

   K = matrix_group(S[1],S[2],S[1]*S[2])
   x = matrix(F,2,2,[2,z,0,2])
   @test x in K
   @test isdefined(K,:X)
   @test isdefined(x,:X)

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
   @test x*y.elm isa fq_nmod_mat
   @test (x*y).elm==x.elm*y
   @test G(x*y.elm)==x*y   
   @test matrix(x)==x.elm


   xg = GAP.Globals.Random(G.X)
   yg = GAP.Globals.Random(G.X)
   pg = MatrixGroupElem(G, xg*yg)
   @test pg==MatrixGroupElem(G,G.mat_iso(xg))*MatrixGroupElem(G,G.mat_iso(yg))

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
   @test isdefined(S,:mat_iso)
   @test S.mat_iso==G.mat_iso
   O = GO(1,2,F)
   H = intersect(S,O)[1]
   @test H==SO(1,2,F)
   @test isnormal(O,H)
   @test index(O,H)==2
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
   @test representative(cc)==x
   @test parent(representative(cc))==G
   @test length(cc)==index(G,C)

   cc = conjugacy_class(G,H)
   @test H^G[2] in elements(cc)
   @test representative(cc)==H
   @test length(cc)==index(G,normalizer(G,H)[1])


   G = GL(2,3)
   @test length(conjugacy_classes(G))==8
   @test length(conjugacy_classes_subgroups(G))==16
   @test length(conjugacy_classes_maximal_subgroups(G))==3
end

