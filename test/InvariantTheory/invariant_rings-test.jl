@testset "InvariantRings" begin
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])

  RG = invariant_ring(M1, M2)
  @test coefficient_ring(RG) == K
  @test !ismodular(RG)

  R = Oscar.polynomial_ring(RG)
  @test Oscar.reynolds_operator(RG, gens(R)[3]^3) == gens(R)[3]^3
  @test Oscar.reynolds_operator(RG, gens(R)[1]) == zero(R)
  mol = Oscar.molien_series(RG)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == (-t^6 - t^3 - 1)//(t^12 - 2t^9 + 2t^3 - 1)

  @test length(Oscar.invariant_basis(RG, 1)) == 0
  @test length(Oscar.invariant_basis(RG, 3)) == 3

  primaries = primary_invariants(RG)
  @test dim(ideal(R, primaries)) == 0
  for f in primaries
    @test Oscar.reynolds_operator(RG, f) == f
  end

  secondaries = secondary_invariants(RG)
  for f in secondaries
    @test Oscar.reynolds_operator(RG, f) == f
  end
  irrs = Oscar.irreducible_secondary_invariants(RG)
  for f in irrs
    @test Oscar.reynolds_operator(RG, f) == f
  end

  K = GF(3)
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 2, 0, 0, 0, 0, 2 ])
  M2 = matrix(K, 3, 3, [ 2, 0, 0, 0, 2, 0, 0, 0, 2 ])

  RG = invariant_ring(M1, M2)
  @test coefficient_ring(RG) == K
  @test !ismodular(RG)

  R = Oscar.polynomial_ring(RG)
  @test Oscar.reynolds_operator(RG, gens(R)[3]^2) == gens(R)[3]^2
  @test Oscar.reynolds_operator(RG, gens(R)[1]) == zero(R)
  mol = Oscar.molien_series(RG)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == (-t^4 - 1)//(t^8 - 2t^6 + 2t^2 - 1)

  @test length(Oscar.invariant_basis(RG, 1)) == 0
  @test length(Oscar.invariant_basis(RG, 2)) == 2

  primaries = primary_invariants(RG)
  @test dim(ideal(R, primaries)) == 0
  for f in primaries
    @test Oscar.reynolds_operator(RG, f) == f
  end

  secondaries = secondary_invariants(RG)
  for f in secondaries
    @test Oscar.reynolds_operator(RG, f) == f
  end
  irrs = Oscar.irreducible_secondary_invariants(RG)
  for f in irrs
    @test Oscar.reynolds_operator(RG, f) == f
  end


end
