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

  @test length(basis(RG, 1)) == 0
  @test length(basis(RG, 1, :reynolds)) == 0
  @test length(basis(RG, 1, :linear_algebra)) == 0
  @test length(basis(RG, 3)) == 3
  @test length(basis(RG, 3, :reynolds)) == 3
  @test length(basis(RG, 3, :linear_algebra)) == 3

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

  @test length(basis(RG, 1)) == 0
  @test length(basis(RG, 1, :reynolds)) == 0
  @test length(basis(RG, 1, :linear_algebra)) == 0
  @test length(basis(RG, 2)) == 2
  @test length(basis(RG, 2, :reynolds)) == 2
  @test length(basis(RG, 2, :linear_algebra)) == 2

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

  # S4
  G = matrix_group(matrix(QQ, [-1 1 0 0;
                               -1 0 1 0;
                               -1 0 0 1;
                               -1 0 0 0]),
                   matrix(QQ, [0 1 0 0;
                               1 0 0 0;
                               0 0 1 0;
                               0 0 0 1]))
  I = invariant_ring(G)
  S, t = QQ["t"]
  m = @inferred molien_series(S, I)
  @test m == 1//((1 - t^2)*(1 - t^3)*(1 - t^4)*(1 - t^5))
  @test m == Oscar._molien_series_via_singular(S, I)

  gl = general_linear_group(4, 5)
  gapmats = [GAP.Globals.PermutationMat(elm.X, 4, GAP.Globals.GF(5))
             for elm in gens(symmetric_group(4))]
  s4 = sub(gl, [MatrixGroupElem(gl, preimage(gl.mat_iso, x), x) for x in gapmats])[1]
  I = invariant_ring(s4)
  m = @inferred molien_series(S, I)
  @test m == 1//((1 - t)*(1 - t^2)*(1 - t^3)*(1 - t^4))

  F = GF(3)
  I = invariant_ring(-identity_matrix(F, 2))
  m = @inferred molien_series(S, I)
  @test m == (t^2 + 1)//(t^4 - 2*t^2 + 1)
  @test m == Oscar._molien_series_via_singular(S, I)
end
