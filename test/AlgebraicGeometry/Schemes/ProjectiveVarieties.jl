@testset "Projective Varieties" begin
  P,_ = polynomial_ring(QQ, [:x, :y, :z])
  G,(x,y,z) = grade(P)
  p = x^2+y^2+z^2
  X = projective_variety(p)
  @test_throws ErrorException projective_variety(p^2)
  Y = projective_scheme(G, ideal([p^2]))
  @test X != Y
end
