@testset "Spectra" begin
  An = affine_space( QQ, 4 )
  @test isa(An, Spec)
  @test base_ring(An) == QQ
  @test isa(ambient_ring(An), MPolyRing)
  @test isa(defining_ideal(An), MPolyIdeal)
  @test length(denoms(An)) == 0

  R = ambient_ring(An)
  # test that we are caching.
  @test R == ambient_ring(An)
  x = gens(R)
  f = x[1] + x[2]^2 - x[3]^3
  I = ideal(R, [f, x[1]^2 - 1])
  @test isa(subscheme(An, I), Spec)
  X = subscheme(An, f)
  @test isa(X, Spec)
  @test isa(Spec(X), Spec)
  @test isa(Spec(ambient_ring(X)), Spec)
  @test isa(Spec(defining_ideal(X)), Spec)
  Q, p = quo(R, I)
  @test isa(Spec(Q), Spec)
  println(X)
  set_name!(X, "X")
  println(X)
end
