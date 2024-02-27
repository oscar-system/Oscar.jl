@testset "Projective Varieties" begin
  G, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  p = x^2 + y^2 + z^2
  X = variety(p)
  @test_throws ErrorException variety(p^2)
  Y = proj(G, ideal([p^2]))
  @test X != Y
  Y1 = variety(ideal([p]),check=false)
  Y1 = variety(ideal([p]),check=true)
  @inferred is_irreducible(Y1)
  variety(G)
  Q,_ = quo(G, ideal([p]))
  variety(Q)
end
