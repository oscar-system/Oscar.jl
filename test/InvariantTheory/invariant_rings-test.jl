@testset "InvariantRings" begin
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])
  RG0 = invariant_ring(M1, M2)

  # Explicitely call the other constructors
  invariant_ring([ M1, M2 ])
  invariant_ring(K, [ M1, M2 ])
  invariant_ring(matrix_group([ M1, M2 ]))

  F = GF(3)
  N1 = matrix(F, 3, 3, [ 0, 1, 0, 2, 0, 0, 0, 0, 2 ])
  N2 = matrix(F, 3, 3, [ 2, 0, 0, 0, 2, 0, 0, 0, 2 ])
  RGp = invariant_ring(N1, N2) # char p, non-modular

  N3 = matrix(F, 2, 2, [ 1, 1, 0, 1 ])
  RGm = invariant_ring(N3) # charp, modular

  @test coefficient_ring(RG0) == K
  @test coefficient_ring(RGp) == F
  @test coefficient_ring(RGm) == F

  @test !is_modular(RG0)
  @test !is_modular(RGp)
  @test is_modular(RGm)

  R0 = polynomial_ring(RG0)
  Rp = polynomial_ring(RGp)
  Rm = polynomial_ring(RGm)

  @test reynolds_operator(RG0, gens(R0)[3]^3) == gens(R0)[3]^3
  @test reynolds_operator(RG0, gens(R0)[1]) == zero(R0)

  @test reynolds_operator(RGp, gens(Rp)[3]^2) == gens(Rp)[3]^2
  @test reynolds_operator(RGp, gens(Rp)[1]) == zero(Rp)

  @test_throws AssertionError reynolds_operator(RGm, gens(Rm)[1])

  @test length(basis(RG0, 1)) == 0
  @test length(basis(RG0, 1, :reynolds)) == 0
  @test length(basis(RG0, 1, :linear_algebra)) == 0
  @test length(basis(RG0, 3)) == 3
  @test length(basis(RG0, 3, :reynolds)) == 3
  @test length(basis(RG0, 3, :linear_algebra)) == 3

  @test length(basis(RGp, 1)) == 0
  @test length(basis(RGp, 1, :reynolds)) == 0
  @test length(basis(RGp, 1, :linear_algebra)) == 0
  @test length(basis(RGp, 2)) == 2
  @test length(basis(RGp, 2, :reynolds)) == 2
  @test length(basis(RGp, 2, :linear_algebra)) == 2

  @test length(basis(RGm, 1)) == 1
  @test length(basis(RGm, 1, :linear_algebra)) == 1
  @test_throws AssertionError basis(RGm, 1, :reynolds)

  mol = molien_series(RG0)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == (-t^6 - t^3 - 1)//(t^12 - 2t^9 + 2t^3 - 1)

  mol = molien_series(RGp)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == (-t^4 - 1)//(t^8 - 2t^6 + 2t^2 - 1)

  fund_invars = fundamental_invariants(RG0)
  for f in fund_invars
    @test reynolds_operator(RG0, f) == f
  end
  fund_invars2 = Oscar.fundamental_invariants_via_minimal_subalgebra(RG0)
  for f in fund_invars2
    @test reynolds_operator(RG0, f) == f
  end

  fund_invars = fundamental_invariants(RGp)
  for f in fund_invars
    @test reynolds_operator(RGp, f) == f
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

  gl = general_linear_group(4, 5)
  gapmats = [GAP.Globals.PermutationMat(elm.X, 4, GAP.Globals.GF(5))
             for elm in gens(symmetric_group(4))]
  s4 = sub(gl, [MatrixGroupElem(gl, x) for x in gapmats])[1]
  I = invariant_ring(s4)
  m = @inferred molien_series(S, I)
  @test m == 1//((1 - t)*(1 - t^2)*(1 - t^3)*(1 - t^4))

  F = GF(3)
  I = invariant_ring(-identity_matrix(F, 2))
  m = @inferred molien_series(S, I)
  @test m == (t^2 + 1)//(t^4 - 2*t^2 + 1)
end
