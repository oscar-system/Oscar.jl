@testset "reflected hypercomplex" begin
  # First create some hyper complex
  S, x = graded_polynomial_ring(QQ, [:x, :y, :z])

  r = length(x)
  d = 5
  F = (is_graded(S) ? graded_free_module(S, degree.(x)) : FreeMod(S, r))
  v = sum(u^d*e for (u, e) in zip(x, gens(F)); init = zero(F))

  K = koszul_complex(Oscar.KoszulComplex, v)

  KK = tensor_product(K, K)

  ref = Oscar.ReflectedComplex(KK)

  @test KK[(2, 3)] === ref[(-2, -3)]
  @test map(KK, 1, (2, 3)) === map(ref, 1, (-2, -3))
  @test Oscar.direction(ref, 1) == :cochain
  @test Oscar.original_complex(ref) === KK
end

