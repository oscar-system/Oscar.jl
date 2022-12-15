@testset "lengths of modules" begin
  R, (x,y,z) = QQ["x", "y", "z"]

  f = x^2 + y^2 - z^2 + 1

  P = ideal(R, f)

  L, _ = localization(R, complement_of_ideal(P))

  F = FreeMod(L, 1)

  I = ideal(L, f^2)*ideal(L, [x,y])^5

  M, _ = quo(F, (I*F)[1])

  @test length(M) == 2
  comp = composition_series(M)
  @test all(x->(parent(x) === M), comp)

  N, _ = quo(F, (ideal(L, [x,y,z])*F)[1])
  @test length(N) == 0

  @test_throws ErrorException length(F)

  A, _ = quo(R, ideal(R, [x^2, y, z]))
  W, _ = localization(A, complement_of_ideal(ideal(R, [x,y,z])))
  F = FreeMod(W, 2)
  @test length(F) == 4 
end

