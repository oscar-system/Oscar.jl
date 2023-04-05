@testset "Projective Algebraic Set" begin
  P,_ = polynomial_ring(QQ, [:x, :y, :z])
  G,(x,y,z) = grade(P)

  A = vanishing_locus(x^2 - y^2)
  @test is_reduced(A)
  @test !is_integral(A)
  @test ideal(A) == ideal(G,[x^2-y^2])

  @test length(irreducible_components(A))==2
end
