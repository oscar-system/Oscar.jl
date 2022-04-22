@testset "Posur" begin
  R, (x,y) = QQ["x", "y"]
  U = MPolyPowersOfElement(x+y)
  S = Localization(U)
  F = FreeMod(S, 2)
  Fb = base_ring_module(F)
  A = S[x//(x+y); y//(x+y)^2]
  B, D = clear_denominators(A)
  @test D*A == B

  #b = S[(x+y)*x y]
  b = MatrixSpace(S, 1, 1)((x+y)*x + 5*y//(x+y)^10)
  (success, v) = has_solution(A, b)
  @test success
  @test v*A == b

  V = MPolyComplementOfPrimeIdeal(ideal(R, [x,y]))
  S = Localization(V)
  A = S[x//(x+y+1); y*(x-5)^3]
  B, D = clear_denominators(A)
  @test D*A == B
  b = MatrixSpace(S, 1, 1)((x+y)*x + 5*y)
  (success, v) = has_solution(A, b)
  @test success
  @test v*A == b
end
