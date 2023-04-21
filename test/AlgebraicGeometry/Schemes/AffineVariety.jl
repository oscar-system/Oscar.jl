@testset "affine variety" begin
  P, (x,y) = polynomial_ring(QQ, [:x,:y])
  X = affine_variety(x, check=false)
  @test is_geometrically_integral(X)
end
