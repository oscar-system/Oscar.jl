@testset "Presentation as affine algebra" begin
  K, a = cyclotomic_field(3, "a")
  M1 = matrix(K, 3, 3, [0, 1, 0, 0, 0, 1, 1, 0, 0])
  M2 = matrix(K, 3, 3, [1, 0, 0, 0, a, 0, 0, 0, -a - 1])

  RG = invariant_ring(M1, M2)
  A, AtoR = affine_algebra(RG)
  @test [AtoR(x) for x in gens(A)] == fundamental_invariants(RG)
  @test is_injective(AtoR)

  RG = invariant_ring(M1, M2)
  A, AtoR = affine_algebra(RG; algo_rels=:linear_algebra)
  @test [AtoR(x) for x in gens(A)] == fundamental_invariants(RG)
  @test is_injective(AtoR)

  # [KS99, Example 17.7]
  K, a = cyclotomic_field(3, "a")
  M = matrix(K, 2, 2, [a, 0, 0, a])
  for algo in [:groebner_basis, :linear_algebra]
    RG = invariant_ring(M)
    A, AtoR = affine_algebra(RG; algo_rels=algo)

    R = polynomial_ring(RG).R
    x = gens(R)

    # Compute the preimages of the primary and irreducible secondary invariants in
    # KS99 under AtoR
    fs = [x[1]^3, x[2]^3, x[1] * x[2]^2, x[1]^2 * x[2]]

    S = base_ring(A)
    y = gens(S)
    Fs = elem_type(S)[]
    for f in fs
      q, r = divrem(f, [g.f for g in fundamental_invariants(RG)])
      @test is_zero(r)
      @test all(g -> total_degree(g) <= 0, q)
      push!(Fs, sum([q[i]([K(1), K(1)]...) * y[i] for i in 1:length(y)]))
    end

    I = ideal(
      S, [Fs[3]^2 - Fs[2] * Fs[4], Fs[4]^2 - Fs[1] * Fs[3], Fs[3] * Fs[4] - Fs[1] * Fs[2]]
    )
    J = modulus(A)
    @test I == J
  end

  # S_4 as matrix group
  M1 = matrix(QQ, 4, 4, [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0])
  M2 = matrix(QQ, 4, 4, [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
  for algo in [:groebner_basis, :linear_algebra]
    RG = invariant_ring(M1, M2)
    A, = affine_algebra(RG; algo_rels=algo)
    @test is_zero(modulus(A))
  end

  # S_4 as permutation group
  for algo in [:groebner_basis, :linear_algebra]
    RG = invariant_ring(symmetric_group(4))
    A, = affine_algebra(RG; algo_rels=algo)
    @test is_zero(modulus(A))
  end

  # Example of an invariant ring which is not Cohen--Macaulay, see Kemper,
  # "Calculating invariant rings of finite groups over arbitrary fields",
  # Example 10.
  G = permutation_group(4, [cperm([1, 2, 3, 4])])
  RG = invariant_ring(GF(2), G)
  A, AtoR = affine_algebra(RG; algo_rels=:groebner_basis)
  @test [AtoR(x) for x in gens(A)] == fundamental_invariants(RG)
  @test is_injective(AtoR)

  RG = invariant_ring(GF(2), G)
  A, AtoR = affine_algebra(RG; algo_rels=:linear_algebra)
  @test [AtoR(x) for x in gens(A)] == fundamental_invariants(RG)
  @test is_injective(AtoR)
end

@testset "Module syzygies" begin
  G = permutation_group(4, [cperm([1, 2, 3, 4])])
  for K in [QQ, GF(2), GF(3)]
    RG = invariant_ring(K, G)
    p_invars = primary_invariants(RG)
    s_invars = secondary_invariants(RG)

    Q, QtoR, StoR = module_syzygies(RG)
    @test ngens(base_ring(Q)) == length(p_invars)
    @test ngens(Q) == length(s_invars)
    S = domain(StoR)
    @test S === base_ring(Q)
    for i in 1:ngens(S)
      @test StoR(gen(S, i)) == p_invars[i]
    end
    for i in 1:rank(Q.F)
      @test QtoR(Q[i]) == s_invars[i]
    end

    if characteristic(K) != 2
      # Non-modular case, so Q must be free
      @test all(is_zero, relations(Q))
    else
      # For GF(2), this is an example of an invariant ring which is not
      # Cohen--Macaulay, see Kemper, "Calculating invariant rings of finite
      # groups over arbitrary fields", Example 10.
      @test !all(is_zero, relations(Q))
      for r in relations(Q)
        @test is_zero(QtoR(r))
      end
    end
  end

  # Another example of a ring which is not Cohen--Macaulay,
  # see DK15, Example 3.6.3.
  for p in [2, 3, 5] # works for an arbitrary prime
    K = GF(p)
    M = matrix(
      K,
      6,
      6,
      [
        1 0 0 0 0 0
        0 1 0 0 0 0
        0 0 1 0 0 0
        1 0 0 1 0 0
        0 1 0 0 1 0
        0 0 1 0 0 1
      ],
    )
    RG = invariant_ring(matrix_group(M))
    p_invars = primary_invariants(RG)
    s_invars = secondary_invariants(RG)

    Q, QtoR, StoR = module_syzygies(RG)
    @test ngens(base_ring(Q)) == length(p_invars)
    @test ngens(Q) == length(s_invars)
    S = domain(StoR)
    @test S === base_ring(Q)
    for i in 1:ngens(S)
      @test StoR(gen(S, i)) == p_invars[i]
    end
    for i in 1:rank(Q.F)
      @test QtoR(Q[i]) == s_invars[i]
    end

    @test !all(is_zero, relations(Q))
    for r in relations(Q)
      @test is_zero(QtoR(r))
    end
  end
end
