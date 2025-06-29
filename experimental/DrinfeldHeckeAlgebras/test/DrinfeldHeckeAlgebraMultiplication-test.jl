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

  @testset "drinfeld-hecke algebra for symplectic doubling of Weyl group B2 over QQ" begin
    # This test was build up using the MAGMA package "CHAMP" by Ulrich Thiel https://github.com/ulthiel/Champ
    # The example is the symplectic reflection algebra for the Weyl group B2
    
    W = matrix_group([matrix(QQ,[-1 2; 0 1]),matrix(QQ,[1 0; 1 -1])])
    G = symplectic_doubling(W)
    A = generic_drinfeld_hecke_algebra(G)
    (t1,t2,t3) = gens(base_ring(A))
    A = evaluate_parameters(A, [-t1,t2,-t3/2])
    
    R = base_ring(A)
    RV = base_algebra(A)
    
    (t,c1,c2) = gens(R)
    (x1,x2,y1,y2) = gens(RV)
    w1 = G[1]
    w2 = G[2]
    
    C = conjugacy_class(G, w2)
    cw2 = collect(C)[2]
    
    @test A[2]*A[1] == A(x1*x2)
    @test A[1]*A[2] == A(x1*x2)
    @test A[3]*A[1] == A(x1 * y1 + t) + c1 * A(w1) + c2 * A(cw2)
    @test A[1]*A[3] == A(x1 * y1)
    @test A[4]*A[1] == A(x1 * y2) + -c2/2 * A(w2) + c2/2 * A(cw2)
    @test A[1]*A[4] == A(x1 * y2)
    @test A[5]*A[1] == -x1 * A(w1)
    @test A[1]*A[5] == x1 * A(w1)
    @test A[6]*A[1] == (x1 + x2) * A(w2)
    @test A[1]*A[6] == x1 * A(w2)
  end
end
