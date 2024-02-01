@testset "views of complexes" begin
  R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
  A = R[x y z; y z x-1]
  I = ideal(R, minors(A, 2))
  Q, _ = quo(R, I)
  L, _ = localization(R, powers_of_element(x-1))
  L1 = FreeMod(L, 1)
  IL1, inc = L(I)*L1
  M = cokernel(inc)
  res, _ = free_resolution(Oscar.SimpleFreeResolution, M);

  K = hom(res, res)
  KV = K[-1, 1:2]

  @test dim(KV) == 1
  @test can_compute_index(KV, 1)
  @test can_compute_index(KV, 2)
  @test !can_compute_index(KV, 3)
  @test !can_compute_index(KV, 0)
  @test Oscar.can_compute_map(KV, 1, (2,))
  @test map(KV, 1, (2,)) === map(K, 2, (-1, 2))

  KV = K[-2, 1:2]

  @test dim(KV) == 1
  @test can_compute_index(KV, 1)
  @test can_compute_index(KV, 2)
  @test !can_compute_index(KV, 3)
  @test !can_compute_index(KV, 0)
  @test Oscar.can_compute_map(KV, 1, (2,))
  @test map(KV, 1, (2,)) === map(K, 2, (-2, 2))

  KV = K[-2:-1, 1]
  @test dim(KV) == 1
  @test can_compute_index(KV, -1)
  @test can_compute_index(KV, -2)
  @test !can_compute_index(KV, -3)
  @test !can_compute_index(KV, 0)
  @test Oscar.can_compute_map(KV, 1, (-1,))
  @test map(KV, 1, (-1,)) === map(K, 1, (-1, 1))

  KV = K[-2:-1, 2]
  @test dim(KV) == 1
  @test can_compute_index(KV, -1)
  @test can_compute_index(KV, -2)
  @test !can_compute_index(KV, -3)
  @test !can_compute_index(KV, 0)
  @test Oscar.can_compute_map(KV, 1, (-1,))
  @test map(KV, 1, (-1,)) === map(K, 1, (-1, 2))

  KV = K[0:0, 1:2]
  @test dim(KV) == 2
  @test can_compute_index(KV, 0, 1)
  @test can_compute_index(KV, 0, 2)
  @test !can_compute_index(KV, 0, 3)
  @test !can_compute_index(KV, 0, 0)
  @test Oscar.can_compute_map(KV, 2, (0, 2))
  @test map(KV, 2, (0, 2)) === map(K, 2, (0, 2))

  KV = K[-10:0, 0:10]
  KV = Oscar.extend_dimension(KV, 3);

  @test dim(KV) == 3
  @test can_compute_index(KV, 0, 1, 0)
  @test can_compute_index(KV, 0, 2, 0)
  @test !can_compute_index(KV, 0, 11, 0)
  @test !can_compute_index(KV, 0, -1, 0)
  @test Oscar.can_compute_map(KV, 2, (0, 2, 0))
  @test map(KV, 2, (0, 2, 0)) === map(K, 2, (0, 2))
end
