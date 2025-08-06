@testset "shifted complexes" begin
  n = 4
  d1 = 3
  d2 = 5
  R, x = polynomial_ring(QQ, n)
  R1 = FreeMod(R, n)
  v = sum(x^d1*e for (x, e) in zip(gens(R), gens(R1)))
  K = Oscar.SimpleComplexWrapper(koszul_complex(v))
  K2 = shift(K, 2) # Creates the same complex, but shifted by 2.

  @test K2[-2] === K[0]
  @test K2[-1] === K[1]
  @test upper_bound(K2) == upper_bound(K) - 2
  @test lower_bound(K2) == lower_bound(K) - 2
  @test all(k->has_index(K2, k) == has_index(K, k+2), -10:10)

  KK = hom(K, K) # A 2-dimensional hypercomplex, actually a double complex
  @test dim(KK) == 2
  KK11 = shift(KK, 1, 1) # The same complex, but shifted by (1, 1)
  @test KK11[0, 0] === KK[1, 1]
  @test is_complete(KK11) == is_complete(KK)
  @test dim(KK11) == 2
  @test all(k->Oscar.direction(KK, k) == Oscar.direction(KK11, k), 1:dim(KK))
end
