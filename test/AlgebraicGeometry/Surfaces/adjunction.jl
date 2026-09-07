@testset "adjunction theory for surfaces" begin
  X = rational_d7_pi4()
  L = adjunction_process(X)
  @test L[1][2][1] == L[1][2][2] == 3
  @test L[1][2][3] == sectional_genus(L[4]) == 1
  X = enriques_d9_pi6()
  L = adjunction_process(X)
  @test L[1][2][4] == 1
end

@testset "fast canonical presentation for adjunction" begin
  X = bordiga()

  Dfast = Oscar._canonical_presentation_matrix_for_adjunction(X)
  Ddirect = Oscar._canonical_presentation_matrix_for_adjunction(X; algorithm = :direct)
  Dmodule = Oscar._canonical_presentation_matrix_for_adjunction(X; algorithm = :module)

  @test Dfast == Ddirect
  @test ncols(Dfast) == ncols(Dmodule)
  @test Oscar._matrix_has_no_quadratic_or_higher_entries(Dfast)
  @test Oscar._matrix_has_no_quadratic_or_higher_entries(Dmodule)

  L = adjunction_process(X, 1; canonical_algorithm = :direct)
  @test L[1][2] == (ZZ(2), ZZ(1), ZZ(0), ZZ(10))
end

@testset "linkage linear strand canonical presentation" begin
  X = rational_d10_pi8()

  D = Oscar._canonical_presentation_matrix_for_adjunction(X; algorithm = :linkage_linear_strand)

  @test base_ring(D) === ambient_coordinate_ring(X)
  @test Oscar._matrix_has_no_quadratic_or_higher_entries(D)
  @test ncols(D) > 0
end

@testset "direct canonical presentation may be nonlinear" begin
  R, x = graded_polynomial_ring(QQ, :x => (1:4))
  X = variety(ideal(R, [x[1]^3 + x[2]^3 + x[3]^3 + x[4]^3]))
  D = Oscar._canonical_presentation_matrix_for_adjunction(X; algorithm = :fast)

  @test !Oscar._matrix_has_no_quadratic_or_higher_entries(D)
  @test D == Oscar._canonical_presentation_matrix_for_adjunction(X; algorithm = :direct)
end

@testset "parametrization" begin
  X = bordiga()
  H = parametrization(X)
  R = domain(H)
  @test total_degree(H(gens(R)[1])) == 4
end
