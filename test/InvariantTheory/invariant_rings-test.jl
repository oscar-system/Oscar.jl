@testset "Invariant rings (for matrix groups)" begin
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])
  RG0 = invariant_ring(M1, M2)

  # Explicitly call the other constructors
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

  @test reynolds_operator(RG0, gens(R0)[3]^3) == reynolds_operator(RG0, gens(R0)[3]^3, trivial_character(group(RG0)))
  @test reynolds_operator(RG0, gens(R0)[1]) == reynolds_operator(RG0, gens(R0)[1], trivial_character(group(RG0)))

  @test reynolds_operator(RGp, gens(Rp)[3]^2) == gens(Rp)[3]^2
  @test reynolds_operator(RGp, gens(Rp)[1]) == zero(Rp)

  @test_throws AssertionError reynolds_operator(RGm, gens(Rm)[1])

  @test length(basis(RG0, 1)) == 0
  @test length(basis(RG0, 1, :reynolds)) == 0
  @test length(basis(RG0, 1, :linear_algebra)) == 0
  @test length(basis(RG0, 1, trivial_character(group(RG0)))) == 0
  @test length(basis(RG0, 3)) == 3
  @test length(basis(RG0, 3, :reynolds)) == 3
  @test length(basis(RG0, 3, :linear_algebra)) == 3
  @test length(basis(RG0, 3, trivial_character(group(RG0)))) == 3

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
  @test molien_series(base_ring(F), RG0, trivial_character(group(RG0))) == mol

  mol = molien_series(RGp)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == (-t^4 - 1)//(t^8 - 2t^6 + 2t^2 - 1)

  # S5 (deleted permutation module)
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

  # S4 (natural permutation module in characteristic 5)
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

@testset "Invariant rings (for permutation groups)" begin
  G = symmetric_group(3)
  RGQ = invariant_ring(G)   # char. 0, over QQ

  K, a = CyclotomicField(3, "a")
  RGK = invariant_ring(K, G)   # char. 0, over K

  F5 = GF(5)
  RGF = invariant_ring(F5, G) # char p, non-modular

  F3 = GF(3)
  RGM = invariant_ring(F3, G)  # char. p, modular

  @test coefficient_ring(RGQ) == QQ
  @test coefficient_ring(RGK) == K
  @test coefficient_ring(RGF) == F5
  @test coefficient_ring(RGM) == F3

  @test !is_modular(RGQ)
  @test !is_modular(RGK)
  @test !is_modular(RGF)
  @test is_modular(RGM)

  RQ = polynomial_ring(RGQ)
  RK = polynomial_ring(RGK)
  RF = polynomial_ring(RGF)
  RM = polynomial_ring(RGM)

  @test reynolds_operator(RGK, gen(RK, 1)^2) == sum(x -> x^2, gens(RK))//3
  @test reynolds_operator(RGK, gen(RK, 1) - gen(RK, 2)) == zero(RK)

  @test reynolds_operator(RGF, gen(RF, 1)^2) == sum(x -> x^2, gens(RF))//3
  @test reynolds_operator(RGF, gen(RF, 1) - gen(RF, 2)) == zero(RF)

  @test_throws AssertionError reynolds_operator(RGM, gen(RM, 1))

  @test length(basis(RGK, 1)) == 1
  @test length(basis(RGK, 1, :reynolds)) == 1
  @test length(basis(RGK, 1, :linear_algebra)) == 1
  @test length(basis(RGK, 3)) == 3
  @test length(basis(RGK, 3, :reynolds)) == 3
  @test length(basis(RGK, 3, :linear_algebra)) == 3

  @test length(basis(RGF, 1)) == 1
  @test length(basis(RGF, 1, :reynolds)) == 1
  @test length(basis(RGF, 1, :linear_algebra)) == 1
  @test length(basis(RGF, 3)) == 3
  @test length(basis(RGF, 3, :reynolds)) == 3
  @test length(basis(RGF, 3, :linear_algebra)) == 3

  @test length(basis(RGM, 1)) == 1
  @test length(basis(RGM, 1, :linear_algebra)) == 1
  @test_throws AssertionError basis(RGM, 1, :reynolds)

  mol = molien_series(RGK)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == 1//((1-t^3)*(1-t^2)*(1-t))

  mol = molien_series(RGF)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == 1//((1-t^3)*(1-t^2)*(1-t))

  # S4 (natural permutation module in characteristic 5)
  s4 = symmetric_group(4)
  S, t = QQ["t"]
  I = invariant_ring(GF(5), s4)
  m = @inferred molien_series(S, I)
  @test m == 1//((1 - t)*(1 - t^2)*(1 - t^3)*(1 - t^4))

  S2 = symmetric_group(2)
  RS2 = invariant_ring(S2)
  R = polynomial_ring(RS2)
  x = gens(R)
  F = abelian_closure(QQ)[1]
  chi = Oscar.group_class_function(S2, [ F(sign(representative(c))) for c in conjugacy_classes(S2) ])
  @test reynolds_operator(RS2, x[1] - x[2], chi) == x[1] - x[2]
  @test reynolds_operator(RS2, x[1] + x[2], chi) == zero(R)

  mol = molien_series(RS2)
  F = parent(mol)
  t = gens(base_ring(F))[1]
  @test mol == 1//(t^3 - t^2 - t + 1)
  @test molien_series(base_ring(F), RS2, chi) == t*mol

  @test length(basis(RS2, 1, chi)) == 1
end
