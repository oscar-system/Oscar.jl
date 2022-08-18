@testset "GAlgebra.printing" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = g_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(gens(r)))

  @test length(string(R)) > 2
  @test length(string(x + y)) > 2
end

@testset "GAlgebra.iteration" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = g_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(gens(r)))

  p = -((x*z*y)^6 - 1)

  s = zero(R)
  for t in terms(p)
    s += t
  end
  @test s == p

  s = zero(R)
  for (c, m) in zip(coefficients(p), monomials(p))
    s += c*m
  end
  @test s == p

  s = build_ctx(R)
  for (c, e) in zip(coefficients(p), exponent_vectors(p))
    push_term!(s, c, e)
  end
  @test finish(s) == p
end
