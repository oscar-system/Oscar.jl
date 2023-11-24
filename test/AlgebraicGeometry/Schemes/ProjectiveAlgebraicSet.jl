@testset "Projective Algebraic Set" begin
  G, (x,y,z) = graded_polynomial_ring(QQ, [:x, :y, :z])

  A = algebraic_set(x^2 - y^2)
  @test is_reduced(A)
  @test !is_integral(A)
  @test fat_ideal(A) == ideal(G, [x^2 - y^2])

  @test length(irreducible_components(A)) == 2
  E = algebraic_set(ideal(G, [x, y^2, z]))
  @test is_empty(E)
  @test !is_integral(E)
  @test_throws ErrorException ProjectiveVariety(E)
end
