@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeAlgebraMultiplication" begin
  @testset "drinfeld-hecke algebra for C2 over QQ" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    R = base_ring(A)
    RV = base_algebra(A)
    
    (s,t) = gens(R)
    (x,y) = gens(RV)
    g = G[1]
    
    RVG = A.form.group_algebra
    
    @testset "multiply with scalar" begin
      @test 2 * A(x) == A(2*x)
      @test A(x) * 2 == A(2*x)
    end
    
    @testset "multiply x * y" begin
      @test A(x) * A(y) == A(x*y)
    end
  
    @testset "multiply y * x" begin
      @test A(y) * A(x) == A(x*y - s) - A(t) * g
    end
  
    @testset "multiply x * g" begin
      @test A(x) * A(g) == A(RVG(x) * RVG(g))
    end
  
    @testset "multiply h * x" begin
      @test A(g) * A(x) == A(RVG(-x) * RVG(g))
    end
  
    @testset "multiply xy * xy" begin
      @test A(x*y) * A(x*y) == A(x^2*y^2 - s * x*y) + A(t*x*y) * g
    end
  
    @testset "multiply xg * yg" begin
      @test A(x) * g * A(y) * g == -A(x*y)
    end
  
    @testset "multiply yg * xg" begin
      @test A(y) * g * A(x) * g == -A(x*y - s) + A(t) * g
    end
  end
end
