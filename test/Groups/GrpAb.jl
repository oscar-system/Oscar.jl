@testset "Groups/GrpAb.jl" begin
  A = abelian_group([4])
  f = hom(A, A, [2*gen(A, 1)])
  g = Oscar.restrict_codomain(f)
  @test is_surjective(g)
  @test order(codomain(g)) == 2
  
  A = abelian_group([0, 2])
  f = hom(A, A, [gen(A, 2), zero(A)])
  g = Oscar.restrict_codomain(f)
  @test is_surjective(g)
  @test order(codomain(g)) == 2
end

@testset "describe for FinGenAbGroup" begin
  @test describe(abelian_group(FinGenAbGroup, Int[])) == "0"
  @test describe(abelian_group(FinGenAbGroup, Int[0])) == "Z"
  @test describe(abelian_group(FinGenAbGroup, Int[0, 0])) == "Z^2"
  @test describe(abelian_group(FinGenAbGroup, Int[2])) == "Z/2"
  @test describe(abelian_group(FinGenAbGroup, Int[0, 2])) == "Z/2 + Z"
  @test describe(abelian_group(FinGenAbGroup, Int[0, 0, 2])) == "Z/2 + Z^2"
  @test describe(abelian_group(FinGenAbGroup, Int[2, 4])) == "Z/2 + Z/4"
  @test describe(abelian_group(FinGenAbGroup, Int[0, 2, 4])) == "Z/2 + Z/4 + Z"
  @test describe(abelian_group(FinGenAbGroup, Int[0, 0, 2, 4])) == "Z/2 + Z/4 + Z^2"
end

@testset "group functions for finite FinGenAbGroup" begin
  @testset for para in [ [2, 3, 4], Int[], [2, 4] ]
    G1 = abelian_group(FinGenAbGroup, para)
    iso = isomorphism(PcGroup, G1)
    G2 = codomain(iso)
    primes = [p for (p, e) in factor(order(G1))]

    # elements
    @test [order(iso(x)) for x in gens(G1)] == [order(x) for x in gens(G1)]
    for x in gens(G1), y in gens(G1)
      @test is_one(comm(x, y))
    end

    # properties
    @test is_abelian(G1) == is_abelian(G2)
    @test is_elementary_abelian(G1) == is_elementary_abelian(G2)
    @test is_finite(G1) == is_finite(G2)
    @test is_finitely_generated(G1) == is_finitely_generated(G2)
    @test is_perfect(G1) == is_perfect(G2)
    @test is_pgroup(G1) == is_pgroup(G2)
    @test is_quasisimple(G1) == is_quasisimple(G2)
    @test is_simple(G1) == is_simple(G2)
    @test is_solvable(G1) == is_solvable(G2)
    @test is_sporadic_simple(G1) == is_sporadic_simple(G2)

    # attributes
    @test images(iso, fitting_subgroup(G1)[1]) == fitting_subgroup(G2)
    @test images(iso, frattini_subgroup(G1)[1]) == frattini_subgroup(G2)
    @test is_pgroup_with_prime(G1) == is_pgroup_with_prime(G2)
    @test nilpotency_class(G1) == nilpotency_class(G2)
    @test number_of_conjugacy_classes(G1) == order(G1)
    @test number_of_conjugacy_classes(Int, G1) isa Int
    @test order(G1) == order(G2)
    @test images(iso, socle(G1)[1]) == socle(G2)
    @test images(iso, solvable_radical(G1)[1]) == solvable_radical(G2)
    if is_pgroup(G1) && order(G1) != 1
      @test prime_of_pgroup(G1) == prime_of_pgroup(G2)
    end
    sg1 = small_generating_set(G1)
    @test order(sub(G1, sg1)[1]) == order(G1)
    @test length(sg1) <= length(para)
    @test abelian_invariants(G1) == abelian_invariants(G2)
    @test abelian_invariants(Int, G1) isa Vector{Int}

    # conjugacy classes of elements
    cc = conjugacy_classes(G1)
    @test length(cc) == order(G1)
    @test all(C -> length(C) == 1, cc)
    @test rand(cc[1]) == representative(cc[1])
    @test acting_group(cc[1]) == G1
    for x in gens(G1), y in gens(G1)
      @test is_conjugate(G1, x, y) == (x == y)
      @test is_conjugate_with_data(G1, x, y)[1] == (x == y)
    end
    C = cc[1]
    for H in C
      @test H == representative(C)
    end

    # conjugacy classes of subgroups
    CC = subgroup_classes(G1)
    @test length(CC) == length(subgroup_classes(G2))
    @test length(CC) == length(collect(subgroups(G1)))
    @test all(C -> length(C) == 1, CC)
    @test rand(CC[1]) == representative(CC[1])
    @test acting_group(CC[1]) == G1
    for C1 in CC, x in gens(G1)
      H = representative(C1)
      @test conjugate_group(H, x) == H
      @test H^x == H
    end
    for C1 in CC, C2 in CC
      H = representative(C1)
      K = representative(C2)
      @test is_conjugate(G1, H, K) == (H == K)
      @test is_conjugate_with_data(G1, H, K)[1] == (H == K)
      @test is_conjugate_subgroup(G1, H, K) == is_subset(K, H)
      @test is_conjugate_subgroup_with_data(G1, H, K) == (is_subset(K, H), zero(G1))
    end
    C = CC[1]
    for H in C
      @test H == representative(C)
    end
    S1 = map(representative, subgroup_classes(G1))
    S2 = map(representative, subgroup_classes(G2))
    @test sort!([order(x) for x in S1]) == sort!([order(x) for x in S2])
    for n in 2:4
      S1 = subgroup_classes(G1, order = n)
      S2 = subgroup_classes(G2, order = n)
      @test length(S1) == length(S2)
    end
    for n in 1:4
      S1 = low_index_subgroup_classes(G1, n)
      S2 = low_index_subgroup_classes(G2, n)
      @test length(S1) == length(S2)
      @test length(S1) == length(collect(low_index_subgroups(G1, n)))
    end
    S1 = maximal_subgroup_classes(G1)
    S2 = maximal_subgroup_classes(G2)
    @test sort!([length(x) for x in S1]) == sort!([length(x) for x in S2])
    @test length(S1) == length(collect(maximal_subgroups(G1)))

    # operations
    x = representative(rand(cc))
    H = representative(rand(CC))
    @test images(iso, core(G1, H)[1]) == core(G2, images(iso, H)[1])
    @test images(iso, normalizer(G1, x)[1]) == normalizer(G2, iso(x))
    @test images(iso, normalizer(G1, H)[1]) == normalizer(G2, images(iso, H)[1])
    @test images(iso, normal_closure(G1, H)[1]) == normal_closure(G2, images(iso, H)[1])

    # operations depending on primes
    for p in primes
      @test images(iso, pcore(G1, p)[1]) == pcore(G2, p)
      @test images(iso, sylow_subgroup(G1, p)[1]) == sylow_subgroup(G2, p)
    end

    # operations depending on sets of primes
    for P in subsets(Set(primes))
      @test [images(iso, representative(C))[1] for C in hall_subgroup_classes(G1, collect(P))] ==
            map(representative, hall_subgroup_classes(G2, collect(P)))
      @test [images(iso, C)[1] for C in hall_subgroups(G1, collect(P))] ==
            collect(hall_subgroups(G2, collect(P)))
    end
    @test issetequal(
      [order(images(iso, S)[1]) for S in hall_system(G1)],
      [order(S) for S in hall_system(G2)]
    )
    @test issetequal([images(iso, S)[1] for S in sylow_system(G1)], sylow_system(G2))
  end
end

@testset "conversions between formats of abelian invariants" begin
  @test Oscar.elementary_divisors_of_vector(Int, []) == []
  @test Oscar.elementary_divisors_of_vector(Int, [0, 3, 2]) == [6, 0]
  @test Oscar.abelian_invariants_of_vector(Int, []) == []
  @test Oscar.abelian_invariants_of_vector(Int, [0, 6]) == [0, 2, 3]
  for i in 1:100
    v = rand(-5:30, 10)
    elab = Oscar.elementary_divisors_of_vector(Int, v)
    abinv = Oscar.abelian_invariants_of_vector(Int, v)
    @test Oscar.elementary_divisors_of_vector(Int, abinv) == elab
    @test Oscar.abelian_invariants_of_vector(Int, elab) == abinv
    @test elementary_divisors(abelian_group([abs(x) for x in v])) == elab
  end
end

@testset "abelian_invariants_schur_multiplier for FinGenAbGroup" begin
  for g in all_small_groups(1:50, is_abelian)
    gg = codomain(isomorphism(FinGenAbGroup, g))
    @test abelian_invariants_schur_multiplier(g) == abelian_invariants_schur_multiplier(gg)
    @test abelian_invariants(schur_multiplier(g)) == abelian_invariants_schur_multiplier(g)
  end

  @test schur_multiplier(PcGroup, abelian_group([2, 3, 4])) isa PcGroup
end
