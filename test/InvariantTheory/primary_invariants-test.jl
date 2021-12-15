@testset "Primary invariants" begin
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])
  RG0 = invariant_ring(M1, M2)

  F = GF(3)
  N1 = matrix(F, 3, 3, [ 0, 1, 0, 2, 0, 0, 0, 0, 2 ])
  N2 = matrix(F, 3, 3, [ 2, 0, 0, 0, 2, 0, 0, 0, 2 ])
  RGp = invariant_ring(N1, N2) # char p, non-modular

  N3 = matrix(F, 2, 2, [ 1, 1, 0, 1 ])
  RGm = invariant_ring(N3) # char p, modular

  for RG in [ RG0, RGp ]
    for algo in [ :optimal_hsop, :radical_containment ]
      invars = primary_invariants(RG, algo)
      @test dim(ideal(polynomial_ring(RG), invars)) == 0
      for f in invars
        @test reynolds_operator(RG, f) == f
      end
    end
  end

  for algo in [ :optimal_hsop, :radical_containment ]
    invars = primary_invariants(RGm, algo)
    @test dim(ideal(polynomial_ring(RGm), invars)) == 0
    actionN3 = Oscar.right_action(polynomial_ring(RGm), N3)
    for f in invars
      @test actionN3(f) == f
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
