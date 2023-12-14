@testset "printing of hypercomplexes" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  S1 = graded_free_module(S, [0])
  I, inc = ideal(S, gens(S))*S1
  M = cokernel(inc)
  res, aug = free_resolution(Oscar.SimpleFreeResolution, M)
  str1 = "$(res)"
  str2 = "$(Oscar.underlying_complex(res))"

  t = tensor_product(res, res)
  str3 = "$(t)"
end
