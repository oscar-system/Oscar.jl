@testset "Examples.ModStdNF" begin
  Qx, x = PolynomialRing(QQ, ["x"])
  I = ideal(Qx, [zero(Qx)])
  G, m = Oscar._compute_standard_basis_with_transform(I)
  @test length(G) == 0
  @test size(m) == (1, 0)

  I = Oscar.DerejeGB.example_1()
  @test length(groebner_basis(I)) == 3
end


