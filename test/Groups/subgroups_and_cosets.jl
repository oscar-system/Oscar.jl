@testset "Subgroups" begin
   S = symmetric_group(7)
   Sx = cperm(S, [1, 2, 3, 4, 5, 6, 7])
   Sy = cperm(S, [1, 2, 3])
   Sz = cperm(S, [1, 2])

   for T in [PermGroup, FPGroup]
     iso = isomorphism(T, S)
     G = codomain(iso)
     x = iso(Sx)
     y = iso(Sy)
     z = iso(Sz)

     H, f = sub(G, [x, y])
     K, g = sub(x, y)
     @test H isa Oscar.sub_type(T)
     @test domain(f) == H
     @test codomain(f) == G
     @test [f(x) for x in gens(H)] == gens(H)
     @test (H, f) == (K, g)
     @test is_subset(K, G)
     flag, emb = is_subgroup(K, G)
     @test flag
     @test g == emb
     @test g == embedding(K, G)
     @test K === domain(emb)
     @test G === codomain(emb)
     @test is_normal_subgroup(H, G)
     H, f = sub(G, [x, z])
     @test H == G
     @test f == id_hom(G)
   end

   G = symmetric_group(7)
   x = cperm(G, [1, 2, 3, 4, 5, 6, 7])
   y = cperm(G, [1, 2, 3])
   H, f = sub(G, [x, y])
   @test H == alternating_group(7)
   @test !is_subset(symmetric_group(8), G)
   @test_throws ArgumentError embedding(symmetric_group(8), G)

   H=sub(G,[G([2,3,1]),G([2,1])])[1]
   @test H != symmetric_group(3)
   @test is_isomorphic(H, symmetric_group(3))
   @test Vector(H[1])==[2,3,1,4,5,6,7]
   @test Vector(symmetric_group(3)(H[1]))==[2,3,1]

   G = symmetric_group(4)
   A = alternating_group(4)
   L = collect(subgroups(G))
   @test length(L)==30
   @test L[1] isa PermGroup
   L1 = [x for x in L if is_normal_subgroup(x, G)]
   K = normal_subgroups(G)
   @test length(K)==4
   for H in L1
      @test H in K
   end
   @test length(maximal_subgroup_classes(G)) == 3
   @test sum(map(length, maximal_subgroup_classes(G))) == 8
   @test any(C -> A in C, maximal_subgroup_classes(G))
   @test maximal_normal_subgroups(G)==[A]
   H = sub(G,[G([3,4,1,2]), G([2,1,4,3])])[1]
   @test minimal_normal_subgroups(G)==[H]
   @test length(characteristic_subgroups(G))==4
   @test H in characteristic_subgroups(G)
   @test is_characteristic_subgroup(H, G)

   H1,h1 = sub(G, gens(symmetric_group(3)))
   H2,h2 = sub(G, gens(alternating_group(4)))
   K,f1,f2 = intersect(H1,H2)
   @test domain(f1)==K
   @test domain(f2)==K
   @test codomain(f1)==H1
   @test codomain(f2)==H2
   @test f1*h1==f2*h2
   @test degree(K)==degree(G)
   @test (K,f1*h1)==sub(G, gens(alternating_group(3)))

   @test derived_subgroup(G)[1] == alternating_group(4)
   L = derived_series(G)
   @test L[1] == G
   @test L[2] == alternating_group(4)
   @test L[3] == sub(G, [cperm([1,2],[3,4]), cperm([1,3],[2,4])])[1]
   @test L[4] == sub(G, [one(G)])[1]

   G = symmetric_group(4)
   A = alternating_group(4)
   @test_throws ArgumentError sub(A, gens(G))
end

@testset "Centralizers and Normalizers in Sym(n)" begin
   G = symmetric_group(6)
   x = cperm(G,[1,2])
   y = cperm(G,[1,2,3,4])

   Cx = centralizer(G,x)[1]
   Cy = centralizer(G,y)[1]
   @testset for i in 1:ngens(Cx)
      @test Cx[i]*x == x*Cx[i]
   end
   @testset for i in 1:ngens(Cy)
      @test Cy[i]*y == y*Cy[i]
   end
   notx = setdiff(G,Cx)
   noty = setdiff(G,Cy)
   @testset for i in 1:3
      z=rand(notx)
      @test z*x != x*z
      z=rand(noty)
      @test z*y != y*z
   end
   @test x in Cx
   @test one(G) in Cx
   @test is_isomorphic(Cx, direct_product(symmetric_group(2),symmetric_group(4)))
   @test is_isomorphic(Cy, direct_product(sub(G,[y])[1], symmetric_group(2)))

   Nx = normalizer(G,Cx)[1]
   Ny = normalizer(G,Cy)[1]
   @test is_normal_subgroup(Cx, Nx)
   @test is_normal_subgroup(Cy, Ny)
   notx = setdiff(G,Nx)
   noty = setdiff(G,Ny)
   @testset for i in 1:3
      z=rand(notx)
      @test Cx^z != Cx
      z=rand(noty)
      @test Cy^z != Cy
   end

   CCx = centralizer(G,Cx)[1]
   CCy = centralizer(G,Cy)[1]
   @testset for i in 1:ngens(Cx), j in 1:ngens(CCx)
      @test Cx[i]*CCx[j] == CCx[j]*Cx[i]
   end
   @testset for i in 1:ngens(Cy), j in 1:ngens(CCy)
      @test Cy[i]*CCy[j] == CCy[j]*Cy[i]
   end
   notx = setdiff(G,CCx)
   noty = setdiff(G,CCy)
   @testset for i in 1:3
      z=rand(notx)
      @test 1 in [z*k != k*z for k in gens(Cx)]
      z=rand(noty)
      @test 1 in [z*k != k*z for k in gens(Cy)]
   end

   Q=quaternion_group(16)
   H=sub(Q,[Q[1]])[1]
   C=core(Q,H)[1]
   @test is_normal_subgroup(C, Q)
   @test order(C)==2
   S=symmetric_group(4)
   P2=pcore(S,2)[1]
   @test order(P2)==4
   @test is_normal_subgroup(P2, S)
   P3=pcore(S,3)[1]
   @test order(P3)==1
   @test_throws ArgumentError pcore(S,4)
end

@testset "Cosets" begin
  
   G = dihedral_group(8)
   H, mH = center(G)
  
   @test index(G, H) == 4
  
   C = right_coset(H, G[1])
   @test group(C) == G
   @test acting_group(C) == H
   @test representative(C) in C
   @test all(x -> x/G[1] in H, C)
   @test is_right(C)
   @test order(C) == length(collect(C))
   @test order(C) isa ZZRingElem
   @test order(Int, C) isa Int
   @test AbstractAlgebra.PrettyPrinting.repr_terse(C) == "Right coset of a group"
  
   T = right_transversal(G, H)
   @test length(T) == index(G, H)
   @test AbstractAlgebra.PrettyPrinting.repr_terse(T) == "Right transversal of groups"
   @test_throws ArgumentError right_transversal(H, G)
   @test_throws ArgumentError left_transversal(H, G)

   G, _ = stabilizer(symmetric_group(6), 6)  # smaller than the natural parent
   @testset "set comparison for cosets in PermGroup" begin
      x = G(cperm([1,2,3]))
      y = G(cperm([1,4,5]))
      z = G(cperm([1,2],[3,4]))
      H = sub(G,[y])[1]
      K = sub(G,[z])[1]

      @test_throws ArgumentError rc = right_coset(H,K(z))
      @test_throws ArgumentError rc = left_coset(H,K(z))
      @test_throws ArgumentError dc = double_coset(H,K(z),K)
      @test_throws ArgumentError dc = double_coset(K,K(z),H)
      rc = right_coset(H,x)
      lc = left_coset(H,x)
      dc = double_coset(H,x,K)
      @test rc==H*x
      @test rc == rc * one(H)
      @test lc==x*H
      @test dc==H*x*K
      @test acting_group(rc) == H
      @test acting_group(lc) == H
      @test left_acting_group(dc) == H
      @test right_acting_group(dc) == K
      @test representative(rc) == x
      @test representative(lc) == x
      @test representative(dc) == x
      @test length(collect(rc)) == length(rc)
      @test Set(collect(rc)) == Set(h*x for h in H)
      @test Set(collect(lc)) == Set(x*h for h in H)
      @test Set(h for h in dc) == Set(h*x*k for h in H for k in K)
      @test order(rc) == 3
      @test order(dc) == 6
      r1=rand(rc)
      r2=rand(dc)
      @test r1 in rc
      @test r2 in dc
      @test r1 in dc
      @test issubset(rc,dc)
      @test issubset(left_coset(K,x),dc)
      @test !is_bicoset(rc)

      @inferred GapObj(rc)
      @inferred GapObj(lc)
      @inferred GapObj(dc)

      @test rc == H*x
      @test lc == x*H
      @test dc == H*x*K
      @test rc*y == H*(x*y)
      @test lc*y == (x*y)*H^y
      @test y*rc == H^(y^-1)*(y*x)
      @test y*lc == (y*x)*H
      @test rc*lc == double_coset(H,x^2,H)
   end

   H = sub(G, [cperm(G,[1,2,3]), cperm(G,[2,3,4])])[1]
   L = right_cosets(G,H)
   @test AbstractAlgebra.PrettyPrinting.repr_terse(L) == "Right cosets of groups"
   @test_throws ArgumentError right_cosets(H, G)
   T = right_transversal(G,H)
   @test length(L)==10
   @test length(T)==10
   @test T[end] == T[length(T)]
   @test T == collect(T)
   @testset for t in T
       @test sum([t in l for l in L]) == 1
   end
   @testset for i in 1:length(T)
       @test findfirst(isequal(T[i]), T) == i
   end
   rc = L[1]
   r = representative(rc)
   rc1 = right_coset(H, H[1]*r)
   @test representative(rc1) != representative(rc)
   @test rc1 == rc
   L = left_cosets(G,H)
   @test_throws ArgumentError left_cosets(H, G)
   T = left_transversal(G,H)
   @test length(L)==10
   @test length(T)==10
   @test T[end] == T[length(T)]
   @test T == collect(T)
   @testset for t in T
       @test sum([t in l for l in L]) == 1
   end
   @testset for i in 1:length(T)
       @test findfirst(isequal(T[i]), T) == i
   end
   lc = L[1]
   r = representative(lc)
   lc1 = left_coset(H, r*H[1])
   @test representative(lc1) != representative(lc)
   @test lc1 == lc
   K = sub(G, gens(symmetric_group(3)) )[1]
   x = G([2,3,4,5,1])
   dc = double_coset(H,x,K)
   dc1 = double_coset(H, H[1]*x, K)
   @test group(dc) == G
   @test representative(dc) != representative(dc1)
   @test dc == dc1
   L = double_cosets(G,H,K)
   @test length(L)==3
   @test left_acting_group(L[1])==H
   @test Set([G([4,5,1,2,3]), G([3,4,5,1,2]), G([2,3,4,5,1])])==Set(intersect(dc, sub(G,[x])[1]) )
   lc = left_coset(H,one(G))
   @test Set(intersect(lc,H))==Set(H)
   lc = left_coset(H,x)
   @test intersect(lc,H)==[]
end

@testset "Double cosets" begin
   G = symmetric_group(5)
   x = G(cperm([1, 2, 3]))
   y = G(cperm([1, 4, 5]))
   z = G(cperm([1, 2], [3, 4]))
   H = sub(G, [y])[1]
   K = sub(G, [z])[1]

   dc = double_coset(H, x, K)
   @test !isassigned(dc.X)
   @test !isassigned(dc.size)
   GapObj(dc)
   @test isassigned(dc.X)
   order(dc)
   @test isassigned(dc.size)
   @test GapObj([dc]; recursive = true) isa GapObj
end

@testset "Predicates for groups" begin
   @test !is_simple(alternating_group(4))
   @test is_simple(alternating_group(5))
   @test is_simple(quo(SL(4,3), center(SL(4,3))[1])[1])

   @test is_almost_simple(symmetric_group(5))
   @test !is_simple(symmetric_group(5))

   @test is_perfect(alternating_group(5))
   @test !is_perfect(alternating_group(4))
   
   @test !is_pgroup(alternating_group(4))
   @test is_pgroup(alternating_group(3))
   @test is_pgroup(quaternion_group(8))
   @test is_pgroup(alternating_group(1))

   @test prime_of_pgroup(alternating_group(3)) == 3
   @test prime_of_pgroup(quaternion_group(8)) == 2

   @test is_solvable(alternating_group(4))
   @test !is_solvable(alternating_group(5))
   @test !is_nilpotent(symmetric_group(4))
   @test !is_supersolvable(symmetric_group(4))
   @test is_nilpotent(quaternion_group(8))
   @test is_supersolvable(quaternion_group(8))
   @test nilpotency_class(quaternion_group(8))==2
   @test_throws ArgumentError nilpotency_class(symmetric_group(4))
end

@testset "Schur multiplier" begin
   @test abelian_invariants(schur_multiplier(cyclic_group(1))) == []
   @test abelian_invariants(schur_multiplier(cyclic_group(5))) == []
   @test abelian_invariants(schur_multiplier(symmetric_group(4))) == [2]
   @test abelian_invariants(schur_multiplier(alternating_group(6))) == [2, 3]

   @test schur_multiplier(symmetric_group(4)) isa FinGenAbGroup
   @test schur_multiplier(PcGroup, symmetric_group(4)) isa PcGroup
end

@testset "Schur cover" begin
   @test order(schur_cover(symmetric_group(4))[1]) == 48
   @test order(schur_cover(alternating_group(5))[1]) == 120
   @test order(schur_cover(dihedral_group(12))[1]) == 24

   @test schur_cover(symmetric_group(4))[1] isa FPGroup
   @test schur_cover(PcGroup, symmetric_group(4))[1] isa PcGroup
   @test schur_cover(PermGroup, alternating_group(5))[1] isa PermGroup
end

@testset "Sylow and Hall subgroups" begin
   G = symmetric_group(4)

   P = sylow_subgroup(G,2)[1]
   @test order(P)==8
   @test is_isomorphic(P,dihedral_group(8))
   P = sylow_subgroup(G,3)[1]
   @test order(P)==3
   @test is_conjugate_with_data(G, P, sub(G, [cperm(1:3)])[1])[1]
   P = sylow_subgroup(G,5)[1]
   @test P==sub(G,[one(G)])[1]
   @test_throws ArgumentError P=sylow_subgroup(G,4)

   G = cyclic_group(210)
   g = G[1]
   @testset for i in [2,3,5,7]
      @test sylow_subgroup(G,i) == sub(G,[g^(210 ÷ i)])
   end
   L = [[2],[3],[5],[7],[2,3],[2,5],[2,7],[3,5],[3,7],[5,7],[2,3,5],[2,3,7],[2,5,7],[3,5,7],[2,3,5,7]]
   @testset for l in L
      h = hall_subgroup_classes(G, l)
      @test length(h) == 1
      @test representative(h[1]) == sub(G,[g^(210÷lcm(l))])[1]
   end
   h = hall_subgroup_classes(G, Int64[])
   @test length(h) == 1
   @test representative(h[1]) == sub(G, [one(G)])[1]
   @test length(hall_subgroup_classes(symmetric_group(5), [2, 5])) == 0
   @test_throws ArgumentError hall_subgroup_classes(G, [4])

   L = sylow_system(G)
   Lo = [order(l) for l in L]
   @test length(Lo)==length(factor(order(G)))
   @test prod(Lo) == order(G)
   @test [is_prime(is_perfect_power_with_data(l)[2]) for l in Lo] == [1 for i in 1:length(L)]
   L = complement_system(G)
   Lo = [index(G,l) for l in L]
   @test length(Lo)==length(factor(order(G)))
   @test prod(Lo) == order(G)
   @test [is_prime(is_perfect_power_with_data(l)[2]) for l in Lo] == [1 for i in 1:length(L)]

   L = hall_system(symmetric_group(4))
   @test is_subset(L[1], symmetric_group(4))
   @test Set(order(H) for H in L)==Set(ZZRingElem[1,3,8,24])
   @test_throws ArgumentError hall_system(symmetric_group(5))
   
end

@testset "Complement classes" begin
   # solvable group
   G = symmetric_group(4)
   N = pcore(G, 2)[1]
   C = @inferred complement_classes(G, N)
   @test length(C) == 1

   # nonsolvable factor group
   G = special_linear_group(2, 5)
   N = center(G)[1]
   C = @inferred complement_classes(G, N)
   @test length(C) == 0

   # nonsolvable normal subgroup
   G = symmetric_group(6)
   N = derived_subgroup(G)[1]
   C = @inferred complement_classes(G, N)
   @test length(C) == 2

   # both normal subgroup and factor group nonsolvable:
   # check that GAP throws an error
   # (if not then perhaps a statement in the documentation of
   # `complement_classes` can be changed)
   G = alternating_group(5)
   W = wreath_product(G, G)
   N = kernel(canonical_projection(W))[1]
   @test_throws ErrorException complement_classes(W, N)

   # pc group, with complements
   G = PcGroup(symmetric_group(4))
   N = pcore(G, 2)[1]
   C = @inferred complement_classes(G, N)
   @test length(C) == 1

   # pc group, without complements
   G = dihedral_group(8)
   N = center(G)[1]
   C = @inferred complement_classes(G, N)
   @test length(C) == 0
end

@testset "Some specific subgroups" begin
   G = GL(2,3)
   S = symmetric_group(4)

   @test order(fitting_subgroup(G)[1])==8
   @test fitting_subgroup(S)==sub(S,[S([3,4,1,2]), S([4,3,2,1])])
   @test frattini_subgroup(S)==sub(S,[one(S)])
   @test frattini_subgroup(G)[1]==intersect(collect(maximal_subgroups(G)))[1]
   @test frattini_subgroup(G)==center(G)
   @test is_characteristic_subgroup(center(G)[1], G)
   @test socle(G)==frattini_subgroup(G)
   @test socle(S)==fitting_subgroup(S)   
   @test solvable_radical(S)[1]==S
   S = symmetric_group(5)
   @test solvable_radical(S)==sub(S,[one(S)])
end
