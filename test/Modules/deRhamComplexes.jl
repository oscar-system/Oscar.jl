@testset "de Rham complexes of polynomial rings and their localizations" begin
  R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w], cached=false)
  X = gens(R)

  dX = Oscar.exterior_derivative.(X)
  @test all(w->Oscar.is_kaehler_differential_module(parent(w))[1], dX)
  w = y*dX[1]
  @test Oscar.exterior_derivative(w) == Oscar.exterior_derivative(-x*dX[2])
  W2_alt = Oscar.kaehler_differentials(R, 2, cached=false)
  @test W2_alt === parent(Oscar.exterior_derivative(w, parent=W2_alt))

  W0 = Oscar.kaehler_differentials(R, 0)
  w = x*W0[1]
  @test parent(Oscar.exterior_derivative(w)) === Oscar.kaehler_differentials(R, 1)
  W1_alt = Oscar.kaehler_differentials(R, 1, cached=false)
  @test W1_alt !== Oscar.kaehler_differentials(R, 1)
  @test parent(Oscar.exterior_derivative(w, parent=W1_alt)) === W1_alt

  U = powers_of_element(x)
  L, loc = localization(R, U)
  W = Oscar.de_rham_complex(L)

  d1 = map(W, 1)
  dX = gens(W[1])
  @test d1(L(1//x)*W[1][2]) == -L(1//x^2)*wedge(dX[1], dX[2])

  I = ideal(R, R[4])
  A, pr = quo(R, I)
  WA1 = Oscar.kaehler_differentials(A)
  @test iszero(WA1[4])
  WA = Oscar.de_rham_complex(A)
  @test iszero(WA[4])
  @test !iszero(WA[3])
  d2 = map(WA, 2)
  @test d2(A(x)*WA[2][4]) == WA[3][1]
  @test "$(WA[2])" == "$(Oscar.is_unicode_allowed() ? "Ω" : "\\Omega")^2($A)"

  I = L(I)
  A, pr = quo(L, I)
  WA1 = Oscar.kaehler_differentials(A)
  @test iszero(WA1[4])
  WA = Oscar.de_rham_complex(A)
  @test iszero(WA[4])
  @test !iszero(WA[3])
  d2 = map(WA, 2)
  @test d2(A(1//x)*WA[2][4]) == -A(1//x^2)*WA[3][1]

  # The graded case
  S, _ = grade(R)
  I = ideal(S, S[4])
  A, pr = quo(S, I)
  @test is_graded(A)
  WA1 = Oscar.kaehler_differentials(A)
  @test is_graded(WA1)
  @test iszero(WA1[4])
  WA = Oscar.de_rham_complex(A)
  @test iszero(WA[4])
  @test !iszero(WA[3])
  d2 = map(WA, 2)
  @test d2(A(x)*WA[2][4]) == WA[3][1]
  @test "$(WA[2])" == "$(Oscar.is_unicode_allowed() ? "Ω" : "\\Omega")^2($A)"
end
