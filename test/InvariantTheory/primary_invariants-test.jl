@testset "Primary invariants" begin
  # Kem99, Example 1
  K, a = CyclotomicField(9, "a")
  M = matrix(K, 2, 2, [ a, 0, 0, -a^3 ])
  RG = invariant_ring(M)
  invars = Oscar.primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 2
  @test sort([ total_degree(f.f) for f in invars ]) == [ 6, 9 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0

  # Kem99, Example 2
  K, a = CyclotomicField(9, "a")
  M = matrix(K, 3, 3, [ a, 0, 0, 0, a^2, 0, 0, 0, a^6 ])
  RG = invariant_ring(M)
  invars = Oscar.primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 3
  @test sort([ total_degree(f.f) for f in invars ]) == [ 3, 5, 9 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0

  # Kem99, p. 183: S_3^4
  M1 = diagonal_matrix([ matrix(QQ, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ]) for i = 1:4 ])
  M2 = diagonal_matrix([ matrix(QQ, 3, 3, [ 0, 1, 0, 0, 0, 1, 1, 0, 0 ]) for i = 1:4 ])
  RG = invariant_ring(M1, M2)
  invars = Oscar.primary_invariants_via_optimal_hsop(RG)
  @test length(invars) == 12
  @test sort([ total_degree(f.f) for f in invars ]) == [ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 ]
  @test dim(ideal(polynomial_ring(RG), invars)) == 0
end
