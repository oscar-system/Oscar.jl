@testset "lengths of modules" begin
  R, (x,y,z) = QQ[:x, :y, :z]

  f = x^2 + y^2 - z^2 + 1

  P = ideal(R, f)

  L, _ = localization(R, complement_of_prime_ideal(P))

  F = FreeMod(L, 1)

  I = ideal(L, f^2)*ideal(L, [x,y])^5

  M, _ = quo(F, (I*F)[1])

  @test length(M) == 2
  comp = composition_series(M)
  @test all(x->(parent(x) === M), comp)

  N, _ = quo(F, (ideal(L, [x,y,z])*F)[1])
  @test length(N) == 0

  M1, _ = quo(F, [f*F[1]])
  res = free_resolution(M1; length=PosInf())
  @test rank(res[1]) == 1
  @test rank(res[2]) == 0

  @test_throws ErrorException length(F)

  A, _ = quo(R, ideal(R, [x^2, y, z]))
  W, _ = localization(A, complement_of_prime_ideal(ideal(R, [x,y,z])))
  F = FreeMod(W, 2)
  @test length(F) == 4

  M2, _ = quo(F, [gens(W)[1]*F[1]])
  res = free_resolution(M2; length=2)
  @test length(res) == 2
  @test length(free_resolution(M2; length=0)) == 0
  @test_throws ErrorException free_resolution(M2; length=-1)
  @test_throws ErrorException free_resolution(M2; length=PosInf())
end

@testset "free resolutions over local rings truncated by length" begin
  # `length = n` stops at the n-th free module, `PosInf()` computes everything
  R, (x, y, z) = QQ[:x, :y, :z]
  L, _ = localization(R, complement_of_point_ideal(R, [0, 0, 0]))
  M = subquotient(free_module(L, 1), L[x; y], L[x^2; y^3; z^4])

  fr = free_resolution(M; length=0)
  @test length(fr) == 0
  @test !is_complete(fr)
  @test [rank(fr[i]) for i in 0:4] == [2, 5, 4, 1, 0]
  @test is_complete(fr)
  @test length(free_resolution(M; length=2)) == 2
  @test length(free_resolution(M)) == 3
  @test length(free_resolution_via_kernels(M, 0)) == 0
  @test_throws ErrorException free_resolution(M; length=-1)
  @test_throws ErrorException free_resolution_via_kernels(M, -1)
end

