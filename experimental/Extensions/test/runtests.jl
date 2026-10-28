using HomotopyContinuation

@testset "Extensions" begin
  @testset "HomotopyContinuation" begin
    Qxy, (x, y) = QQ[:x, :y]
    I = ideal([x^2 - x * y + y^3, y^2 - x - 4//5])
    # commented out until we have a correct dim_numerical
    #@test iszero(dim_numerical(I))
    result = Oscar.solve_numerical(I)
    @test [r.multiplicity for r in result.path_results] == [1, 1, 1, 1]

    J = ideal([x^2 - x * y + y^3])
    ws = witness_set(J)
    @test length(ws.R) == 3
  end

  @testset "HomotopyContinuation over QQBar" begin
    K = algebraic_closure(QQ)
    i = sqrt(K(-1))
    Kxy, (x, y) = K[:x, :y]
    result = Oscar.solve_numerical(ideal([x^2 + 1, y - i * x]))
    sols = HomotopyContinuation.solutions(result)
    @test length(sols) == 2
    @test all(s -> isapprox(s[2], im * s[1]; atol=1e-8), sols)

    # rational coefficients over QQBar behave as over QQ
    Qxy, (u, v) = QQ[:x, :y]
    I = ideal([u^2 - u * v + v^3, v^2 - u - 4//5])
    IK = ideal([x^2 - x * y + y^3, y^2 - x - 4//5])
    @test HomotopyContinuation.nsolutions(Oscar.solve_numerical(IK)) ==
          HomotopyContinuation.nsolutions(Oscar.solve_numerical(I))
  end

  # here for future test
  # @testset "" begin
  #   R, a = polynomial_ring(QQ, :a => (1:2, 1:2))
  #   I = ideal(vec(a[:,1]*transpose(a[:,2])-Array(R[1//2 1; 0 1//2])))
  #   Oscar.dim(I)
  #   dim_numerical(I)
  # end
end
