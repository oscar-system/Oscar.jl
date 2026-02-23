@testset "Extensions" begin
  @testset "HomotopyContinuation" begin
    Qxy, (x, y) = QQ[:x, :y]
    @test_throws ErrorException solve_numerical([x^2 - y, x * y])
    @test_throws ErrorException dim_numerical([x^2 - y, x * y])
  end
end
