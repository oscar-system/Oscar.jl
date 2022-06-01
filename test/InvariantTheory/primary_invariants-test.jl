@testset "Primary invariants (for matrix groups)" begin
  # Char 0
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])
  for algo in [ :optimal_hsop, :successive_algo ]
    RG = invariant_ring(M1, M2)
    invars = primary_invariants(RG, algo)
    @test dim(ideal(polynomial_ring(RG), invars)) == 0
    for f in invars
      @test reynolds_operator(RG, f) == f
    end
  end

  # Char p, non-modular
  F3 = GF(3)
  N1 = matrix(F3, 3, 3, [ 0, 1, 0, 2, 0, 0, 0, 0, 2 ])
  N2 = matrix(F3, 3, 3, [ 2, 0, 0, 0, 2, 0, 0, 0, 2 ])
  for algo in [ :optimal_hsop, :successive_algo ]
    RG = invariant_ring(N1, N2)
    invars = primary_invariants(RG, algo)
    @test dim(ideal(polynomial_ring(RG), invars)) == 0
    for f in invars
      @test reynolds_operator(RG, f) == f
    end
  end

  # Char p, modular
  F9, b = FiniteField(3, 2, "b")
  N3 = matrix(F9, [ 1 0 0 0; b + 1 1 0 0; -1 0 1 0; b 0 -1 1 ])
  N4 = matrix(F9, [ 1 0 0 0; 1 1 0 0; 1 0 1 0; b -b b 1 ])
  for algo in [ :optimal_hsop, :successive_algo ]
    RG = invariant_ring(N3, N4)
    invars = primary_invariants(RG, algo)
    @test dim(ideal(polynomial_ring(RG), invars)) == 0
    actionN3 = Oscar.right_action(polynomial_ring(RG), N3)
    actionN4 = Oscar.right_action(polynomial_ring(RG), N4)
    for f in invars
      @test actionN3(f) == f
      @test actionN4(f) == f
    end
  end

  # Kem99, Example 1
  K, a = CyclotomicField(9, "a")
  M = matrix(K, 2, 2, [ a, 0, 0, -a^3 ])
  RG = invariant_ring(M)
  invars = primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 2
  @test sort([ total_degree(f.f) for f in invars ]) == [ 6, 9 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0

  # Kem99, Example 2
  K, a = CyclotomicField(9, "a")
  M = matrix(K, 3, 3, [ a, 0, 0, 0, a^2, 0, 0, 0, a^6 ])
  RG = invariant_ring(M)
  invars = primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 3
  @test sort([ total_degree(f.f) for f in invars ]) == [ 3, 5, 9 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0

  # Kem99, p. 183: S_3^3
  M1 = diagonal_matrix([ matrix(QQ, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ]) for i = 1:3 ])
  M2 = diagonal_matrix([ matrix(QQ, 3, 3, [ 0, 1, 0, 0, 0, 1, 1, 0, 0 ]) for i = 1:3 ])
  RG = invariant_ring(M1, M2)
  invars = primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 9
  @test sort([ total_degree(f.f) for f in invars ]) == [ 1, 1, 1, 2, 2, 2, 3, 3, 3 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0
end

@testset "Primary invariants (for permutation groups)" begin
  # Char 0
  K, a = CyclotomicField(3, "a")
  G = symmetric_group(4)
  for algo in [ :optimal_hsop, :successive_algo ]
    RG = invariant_ring(K, G)
    invars = primary_invariants(RG, algo)
    @test dim(ideal(polynomial_ring(RG), invars)) == 0
    for f in invars
      @test reynolds_operator(RG, f) == f
    end
  end

  # Char p, non-modular
  F5 = GF(5)
  G = symmetric_group(4)
  for algo in [ :optimal_hsop, :successive_algo ]
    RG = invariant_ring(F5, G)
    invars = primary_invariants(RG, algo)
    @test dim(ideal(polynomial_ring(RG), invars)) == 0
    for f in invars
      @test reynolds_operator(RG, f) == f
    end
  end

  # Char p, modular
  F9 = GF(3, 2)
  for algo in [ :optimal_hsop, :successive_algo ]
    RG = invariant_ring(F9, G)
    invars = primary_invariants(RG, algo)
    @test dim(ideal(polynomial_ring(RG), invars)) == 0
    actions = [Oscar.right_action(polynomial_ring(RG), x) for x in gens(G)]
    for f in invars
      @test all(act -> act(f) == f, actions)
    end
  end

  # Kem99, p. 183: S_3^3
  s3 = symmetric_group(3)
  G = PermGroup(direct_product(s3, s3, s3))
  RG = invariant_ring(G)
  invars = primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 9
  @test sort([ total_degree(f.f) for f in invars ]) == [ 1, 1, 1, 2, 2, 2, 3, 3, 3 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0
end
