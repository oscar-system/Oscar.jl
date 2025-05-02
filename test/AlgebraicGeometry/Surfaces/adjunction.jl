@testset "adjunction theory for surfaces" begin
  X = rational_d7_pi4();
  L = adjunction_process(X)
  @test L[1][2][1] == L[1][2][2] == 3
  @test L[1][2][3] == sectional_genus(L[4]) == 1
  X = enriques_d9_pi6();
  L = adjunction_process(X)
  @test L[1][2][4] == 1
end

@testset "parametrization" begin
  X = bordiga()
  H = parametrization(X)
  R = domain(H)
  @test total_degree(H(gens(R)[1])) == 4
end

