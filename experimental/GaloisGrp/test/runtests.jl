@testset "SolveByRadicals" begin
  Qx, x = QQ["x"]
  K, r = solve(x^3+3*x+5)
  @test absolute_degree(K) == 12
  @test length(r) == 3

  Qt, t = rational_function_field(QQ, "t")
  Qtx, x = Qt["x"]
  F, a = function_field(x^6 + 108*t^2 + 108*t + 27)
  subfields(F)
  G, = galois_group(F)
  @test is_isomorphic(G, symmetric_group(3))
end

