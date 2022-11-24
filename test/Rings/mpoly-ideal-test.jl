@testset "MPolyIdeal.vdim" begin
  r, (x, y) = PolynomialRing(QQ, [:x, :y])
  @test vdim(ideal(r, [x^2+y^2])) == -1
  @test vdim(ideal(r, [x^2+y^2, x^2-y^2])) == 4
end
