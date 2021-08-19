@testset "affine schemes" begin
  An = affine_space( QQ, 4 )

  R = ambient_ring(An)
  # test that we are caching.
  @test R == ambient_ring(An)
  x = gens(R)
  I = ideal([x[1] + x[2]^2 - x[3]^3, x[1]^2 - 1])

  A = Spec(R, I)
  f1 = x[1]^2 + x[2]^2 - 1
  f2 = x[3]
  D1 = localize(A, f1)
  D2 = localize(A, f2)
  ambient_ring(D2)
  k = base_ring(An)
  D12 = localize(D1, f2)
  D21 = localize(D2, f1)
  @test base_ring(D1) == k
  @test base_ring(D2) == k
  @test base_ring(D12) == k
  @test root(D1) == A
  @test root(D12) == A
  @test root(D21) == A
  @test D12 != D21
  I = defining_ideal(D12)
end
