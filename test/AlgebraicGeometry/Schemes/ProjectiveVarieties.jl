@testset "Projective Varieties" begin
  P,(x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  G = grade(P)
  p = x^2+y^2+z^2
  X = projective_variety(p)
  Y = projective_scheme(P, ideal([p^2]))
  @test X != Y
end
