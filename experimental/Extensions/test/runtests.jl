using HomotopyContinuation

@testset "Extensions" begin
  @testset "HomotopyContinuation" begin
    Qxy, (x, y) = QQ[:x, :y]
    I = ideal([x^2 - x * y + y^3, y^2 - x - 4//5])
    @test iszero(dim_numerical(I))
    result = solve_numerical(I)
    @test [r.multiplicity for r in result.path_results] == [1, 1, 1, 1]

    J = ideal([x^2 - x * y + y^3])
    ws = witness_set(J)
    @test length(ws.R) == 3
  end

  @testset "" begin
    R, a = polynomial_ring(QQ, :a => (1:2, 1:2))
    I = ideal(vec(a[:,1]*transpose(a[:,2])-Array(R[1//2 1; 0 1//2])))
    Oscar.dim(I)
    dim_numerical(I)
  end

end
