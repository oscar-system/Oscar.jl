@testset "Projective Algebraic Set" begin
  P,_ = polynomial_ring(QQ, [:x, :y, :z])
  G,(x,y,z) = grade(P)

  A = vanishing_locus(x^2 - y^2)
  @test is_reduced(A)
  @test !is_integral(A)
  @test ideal(A) == ideal(G,[x^2-y^2])

  @test length(irreducible_components(A))==2
  E = vanishing_locus(ideal(G, [x,y^2,z]))
  @test is_empty(E)
  @test !is_integral(E)
  @test_throws ErrorException ProjectiveVariety(E)
end
