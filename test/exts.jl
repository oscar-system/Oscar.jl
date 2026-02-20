using HomotopyContinuation

@testset "Extensions" begin
  @testset "HomotopyContinuation" begin
    Qxy, (x, y) = QQ[:x, :y]
    I = ideal([x^2 - x * y + y^3, y^2 - x - 4//5])
    @test iszero(Oscar.numerical_dim(I))
    result = HomotopyContinuation.solve(I)
    @test [r.multiplicity for r in result.path_results] == [1, 1, 1, 1]

    J = ideal([x^2 - x * y + y^3])
    ws = Oscar.witness_set(J)
    @test length(ws.R) == 3
  end
end
