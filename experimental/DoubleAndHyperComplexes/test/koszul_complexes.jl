@testset "koszul complexes" begin
  S, x = graded_polynomial_ring(QQ, [:x, :y, :z])

  r = length(x)
  d = 5
  F = (is_graded(S) ? graded_free_module(S, degree.(x)) : FreeMod(S, r))
  v = sum(u^d*e for (u, e) in zip(x, gens(F)); init = zero(F))

  K = koszul_complex(Oscar.KoszulComplex, v)
  KK = koszul_complex(v)

  for i in 0:4
    @test matrix(map(K, i)) == matrix(map(KK, i))
  end
  
  # the ungraded case
  S, x = polynomial_ring(QQ, [:x, :y, :z])

  r = length(x)
  d = 5
  F = (is_graded(S) ? graded_free_module(S, degree.(x)) : FreeMod(S, r))
  v = sum(u^d*e for (u, e) in zip(x, gens(F)); init = zero(F))

  K = koszul_complex(Oscar.KoszulComplex, [u^d for u in gens(S)]...)
  KK = koszul_complex(v)

  for i in 0:4
    @test matrix(map(K, i)) == matrix(map(KK, i))
  end
end


