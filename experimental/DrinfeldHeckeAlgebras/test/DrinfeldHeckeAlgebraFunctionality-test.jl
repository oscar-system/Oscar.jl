@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeAlgebraFunctionality" begin
  @testset "evaluate parameters" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    B = evaluate_parameters(A, [-8//7, 4])
    R = base_algebra(B)
    (x,y) = gens(R)
    
    @test B(y) * B(x) == B(x*y + 8//7) - 4 * A[3]
  end

  @testset "evaluate parameter not enough values" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)

    @test_throws ArgumentError evaluate_parameters(A, [-8//7])
  end

  @testset "evaluate parameter too many values" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test_throws ArgumentError evaluate_parameters(A, [2,3,4])
  end

  @testset "evaluate parameter wrong element type" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test_throws ArgumentError evaluate_parameters(A, [1, "hi"])
  end

  @testset "evaluate parameter not parametrized" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = drinfeld_hecke_algebra(G)
    
    @test_throws ArgumentError evaluate_parameters(A, [1, 2])
  end
end
