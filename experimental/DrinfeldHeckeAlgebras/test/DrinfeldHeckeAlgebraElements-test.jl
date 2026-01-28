@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeAlgebraElements" begin
  G = matrix_group(matrix(QQ, [-1 0;0 -1]))
  A = generic_drinfeld_hecke_algebra(G)
  
  R = base_ring(A)
  RV = base_algebra(A)
  
  (s,t) = gens(R)
  (x,y) = gens(RV)
  g = G[1]
  
  T = elem_type(typeof(A))
  
  @testset "cast Integer to Drinfeld-Hecke algebra" begin
    @test A(1) isa T
  end

  @testset "cast Rational to Drinfeld-Hecke algebra" begin
    @test A(1//2) isa T
  end

  @testset "cast QQElem to Drinfeld-Hecke algebra" begin
    @test A(QQ(1//2)) isa T
  end

  @testset "cast parameter to Drinfeld-Hecke algebra" begin
    @test A(s) isa T
  end

  @testset "cast base algebra element to Drinfeld-Hecke algebra" begin
    @test A(x) isa T
  end

  @testset "cast parameter * base algebra element to Drinfeld-Hecke algebra" begin
    @test A(s*x) isa T
  end

  @testset "cast group element to Drinfeld-Hecke algebra" begin
    @test A(g) isa T
  end

  @testset "can add mixed types from the left" begin
    @test 1 + A(g) isa T
    @test 1//2 + A(g) isa T
    @test QQ(1//2) + A(g) isa T
    @test s + A(g) isa T
    @test t + A(g) isa T
    @test g + A(g) isa T
  end

  @testset "can add mixed types from the right" begin
    @test A(g) + 1 isa T
    @test A(g) + 1//2 isa T
    @test A(g) + QQ(1//2) isa T
    @test A(g) + s isa T
    @test A(g) + t isa T
    @test A(g) + g isa T
  end

  @testset "can subtract mixed types from the left" begin
    @test 1 - A(g) isa T
    @test 1//2 - A(g) isa T
    @test QQ(1//2) - A(g) isa T
    @test s - A(g) isa T
    @test t - A(g) isa T
    @test g - A(g) isa T
  end

  @testset "can subtract mixed types from the right" begin
    @test A(g) - 1 isa T
    @test A(g) - 1//2 isa T
    @test A(g) - QQ(1//2) isa T
    @test A(g) - s isa T
    @test A(g) - t isa T
    @test A(g) - g isa T
  end

  @testset "can multiply mixed types from the left" begin
    @test 1 * A(g) isa T
    @test 1//2 * A(g) isa T
    @test QQ(1//2) * A(g) isa T
    @test s * A(g) isa T
    @test t * A(g) isa T
    @test g * A(g) isa T
  end

  @testset "can subtract mixed types from the right" begin
    @test A(g) * 1 isa T
    @test A(g) * 1//2 isa T
    @test A(g) * QQ(1//2) isa T
    @test A(g) * s isa T
    @test A(g) * t isa T
    @test A(g) * g isa T
  end
end
