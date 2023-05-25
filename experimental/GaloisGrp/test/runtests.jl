@testset "SolveByRadicals" begin
  Qx, x = QQ["x"]
  K, r = solve(x^3+3*x+5)
  @test absolute_degree(K) == 12
  @test length(r) == 3
end

