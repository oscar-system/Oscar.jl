@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeAlgebra" begin
  MS = matrix_space(QQ,2,2)
  G = matrix_group(matrix(MS([-1 0;0 -1])))
  e = one(G)
  g = gen(G,1)
      
  @testset "create trivial drinfeld-hecke algebra" begin
    A = drinfeld_hecke_algebra(G)
    @test is_zero(A)
  end

  @testset "create drinfeld-hecke algebra from drinfeld-hecke form" begin
    κ = generic_drinfeld_hecke_form(G)
    
    A = drinfeld_hecke_algebra(κ)
    @test !is_zero(A)
  end

  @testset "create parametrized drinfeld-hecke algebra" begin
    A = generic_drinfeld_hecke_algebra(G)
    @test !is_zero(A)
    @test ngens(base_ring(A)) == 2
  end

  @testset "multiplication" begin
    κ = generic_drinfeld_hecke_form(G)
    κ = evaluate_parameters(κ, [-8//7, 4])
    A = drinfeld_hecke_algebra(κ)
    RG = group_algebra(A)
    R = base_algebra(A)
    
    x = R[1]
    y = R[2]
    
    @testset "multiply with scalar" begin
      @test 2 * A(x) == A(2*x)
      @test A(x) * 2 == A(2*x)
    end
    
    @testset "multiply x * y" begin
      @test A(x) * A(y) == A(x*y)
    end
  
    @testset "multiply y * x" begin
      @test A(y) * A(x) == A(x*y + 8//7) - A(4) * g
    end
  
    @testset "multiply x * g" begin
      xg = A(RG(x) * RG(g))
      @test A(x) * A(g) == xg
      @test A(x) * g == xg
#       @test x * A(g) == xg
    end
  
    @testset "multiply h * x" begin
      gx = A(RG(-x) * RG(g))
      @test A(g) * A(x) == gx
#       @test g * A(x) == gx
#       @test A(g) * x == gx
    end
  
    @testset "multiply xy * xy" begin
      @test A(x*y) * A(x*y) == A(x^2*y^2 + 8//7 * x*y) + A(4*x*y) * g
    end
  
    @testset "multiply xg * yg" begin
      @test A(x) * g * A(y) * g == -A(x*y)
    end
  
    @testset "multiply yg * xg" begin
      @test A(y) * g * A(x) * g == -A(x*y + 8//7) + A(4) * g
    end
  end
end
