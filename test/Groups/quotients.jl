@testset "quo for trivial kernel" begin
#  @testset for G in [symmetric_group(4), special_linear_group(2, 3), special_linear_group(2, 4), free_group(1), abelian_group(PcGroup, [2, 3, 4])]
   @testset for G in [free_group(0), symmetric_group(4), special_linear_group(2, 3), special_linear_group(2, 4), free_group(1), abelian_group(SubPcGroup, [2, 3, 4])]
      subgens = elem_type(G)[]
      F, epi = quo(G, subgens)
      if is_finite(G)
        @test order(F) == order(G)
        @test epi(one(G)) == one(F)
      else
        @test ! is_finite(F)
      end
   end
end

@testset "types of factor groups" begin
   # - `quo` without prescribed type:
   #   same type as the input type if the kernel is trivial,
   #   otherwise `PcGroup` if solvable, `PermGroup` if not
   # - `quo` with prescribed type:
   #    return a group of this type if possible
   G = symmetric_group(4)
   N = pcore(G, 2)[1]
   for subgens in [N, gens(N)]
      @test quo(G, subgens)[1] isa PcGroup
      @test quo(PermGroup, G, subgens)[1] isa PermGroup
   end
   N = trivial_subgroup(G)[1]
   for subgens in [N, gens(N)]
      @test quo(G, subgens)[1] isa PermGroup
      @test quo(PcGroup, G, subgens)[1] isa PcGroup
   end
   G = special_linear_group(2, 5)
   N = center(G)[1]
   for subgens in [N, gens(N)]
      @test quo(G, subgens)[1] isa PermGroup
      @test quo(FPGroup, G, subgens)[1] isa FPGroup
      @test_throws ArgumentError quo(PcGroup, G, subgens)
   end
   N = trivial_subgroup(G)[1]
   for subgens in [N, gens(N)]
      @test quo(G, subgens)[1] isa MatrixGroup
      @test quo(PermGroup, G, subgens)[1] isa PermGroup
      @test_throws ArgumentError quo(PcGroup, G, subgens)
   end
   G = automorphism_group(small_group(12, 5))
   N = trivial_subgroup(G)[1]
   for subgens in [N, gens(N)]
      @test quo(G, subgens)[1] isa AutomorphismGroup
      @test quo(PermGroup, G, subgens)[1] isa PermGroup
   end

   # - `maximal_abelian_quotient` without prescribed type:
   #   same type as the input type if abelian,
   #   otherwise `PcGroup` if finite, `FPGroup` if not
   # - `maximal_abelian_quotient` with prescribed type:
   #    return a group of this type if possible
#  for T in [ PermGroup, PcGroup ]
#TODO: support `abelian_group(PcGroup, [2, 3, 4])`, using GAP's PcpGroup
   for T in [ PermGroup, SubPcGroup ]
      G = abelian_group(T, [2, 3, 4])
      @test maximal_abelian_quotient(G)[1] isa T
      @test maximal_abelian_quotient(PermGroup, G)[1] isa PermGroup
   end
   T = FPGroup
   G = abelian_group(T, [2, 3, 4])
   # Note that GAP chooses a representation for the max. abelian quotient,
   # and this group gets wrapped into an Oscar group.
   # Here GAP chooses a pc group, but with noncanonical generators,
   # Hence we get a `SubPcGroup` not a `PcGroup`.
   # (This is not really satisfactory.)
   @test maximal_abelian_quotient(G)[1] isa SubPcGroup
   @test maximal_abelian_quotient(PermGroup, G)[1] isa PermGroup
   @test maximal_abelian_quotient(FinGenAbGroup, G)[1] isa FinGenAbGroup
   G = symmetric_group(4)
   @test maximal_abelian_quotient(G)[1] isa PcGroup
   @test maximal_abelian_quotient(PermGroup, G)[1] isa PermGroup
   G = special_linear_group(2, 3)
   @test maximal_abelian_quotient(G)[1] isa PcGroup
   @test maximal_abelian_quotient(PermGroup, G)[1] isa PermGroup
   G = free_group(1)
   @test maximal_abelian_quotient(G)[1] isa FPGroup
   @test maximal_abelian_quotient(FPGroup, G)[1] isa FPGroup
   @test_throws MethodError quo(PcGroup, G)

   # `abelian_invariants`
   @test abelian_invariants(free_group(2)) == [0, 0]
   @test abelian_invariants(alternating_group(5)) == []
   @test abelian_invariants(small_group(8, 5)) == [2, 2, 2]
end

@testset "Relators" begin
   # free group
   F = free_group(2)
   @test relators(F) == []

   # f.p. group
   x, y = gens(F)
   rels = [x^3, y^2, comm(x, y)]
   G, epi = quo(F, rels)
   @test relators(G) == rels

   # sub-f.p. group
   G = sylow_subgroup(G, 2)[1]
   rels = relators(G)
   @test is_isomorphic(G, quo(parent(rels[1]), rels)[1])

   # pc group
   G = dihedral_group(12)
   rels = relators(G)
   @test is_isomorphic(G, quo(parent(rels[1]), rels)[1])

   # sub-pc group
   G = sylow_subgroup(G, 2)[1]
   rels = relators(G)
   @test is_isomorphic(G, quo(parent(rels[1]), rels)[1])

   # perm. group
   G = symmetric_group(3)
   rels = relators(G)
   @test is_isomorphic(G, quo(parent(rels[1]), rels)[1])
end

@testset "Finitely presented groups" begin
   F = free_group(2)
   x,y = gens(F)
   @test x == F[1]
   @test y == F[2]
   @test Oscar._is_full_fp_group(GapObj(F))
   @test relators(F) == FPGroupElem[]
   S = sub(F, [gen(F, 1)])[1]
   @test ! Oscar._is_full_fp_group(GapObj(S))

   n=5
   G,f = quo(F, [x^2,y^n,(x*y)^2] )           # dihedral group D(2n)
   @test is_finite(G)
   @test order(G) == 2*n
   @test !is_abelian(G)
   @test is_isomorphic(G, dihedral_group(2*n))
   @test !is_injective(f)
   @test is_surjective(f)
   @test exponent(G) == 2*n
   if order(G[2])==n
      @test G[2]^G[1] == G[2]^-1
   else
      @test G[1]^G[2] == G[1]^-1
   end
   
   G,f = quo(F, [x^n,y^n,comm(x,y)])          # group C(n) x C(n)
   @test is_finite(G)
   @test order(G) == n^2
   @test is_abelian(G)
   @test !is_injective(f)
   @test is_surjective(f)
   @test exponent(G) == n
   @test isone(G[1]^n)
   @test Oscar._is_full_fp_group(GapObj(G))
   @test relators(G)==[x^n,y^n,comm(x,y)]
   S = sub(G, [gen(G, 1)])[1]
   @test ! Oscar._is_full_fp_group(GapObj(S))
   @test_throws ArgumentError quo(S, gens(S))

   @test G([1 => 2, 2 => -3]) == G[1]^2 * G[2]^-3
   @test_throws MethodError S([1 => 2])

   S = symmetric_group(4)
   G,f = quo(S, [cperm(S,[1,3,2])])
   @test order(G) == 2
   @test f(S([1,2,4,3]))==G[1]
   @test f(S([2,1,4,3]))==one(G)
end

@testset "matrix groups" begin
   K, a = cyclotomic_field(3, "a");
   S = matrix(K, [0 0 1; 1 0 0; 0 1 0])
   T = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
   H3 = matrix_group(S, T)
   Z = matrix(K, [-a-1 0 0; 0 -a-1 0; 0 0 -a-1])
   C = sub(H3, [H3(Z; check = false)]; check = false)[1]
   @test !has_is_finite(C)
   Q, pQ = quo(H3, C);
   @test has_is_finite(C)
   @test C == center(H3)[1]
end
