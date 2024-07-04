@testset "printing of hypercomplexes" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  S1 = graded_free_module(S, [0])
  I, inc = ideal(S, gens(S))*S1
  M = cokernel(inc)
  res, aug = free_resolution(Oscar.SimpleFreeResolution, M)
  str1 = "$(res)"
  str2 = "$(Oscar.underlying_complex(res))"
  str4 = "$(aug)"

  strand, inc = Oscar.linear_strand(res)
  str = "$(strand)"
  str = "$(inc)"
  q, pr = cokernel(inc)
  str = "$(q)"
  str = "$(pr)"

  t = tensor_product(res, res)
  str3 = "$(t)"

  A, _ = quo(S, ideal(S, z))
  A1 = FreeMod(A, 1)
  I, inc = ideal(A, [x, y])*A1
  M = cokernel(inc)
  res, aug = free_resolution(Oscar.SimpleFreeResolution, M)
  str = "$(res)"
  str = "$(aug)"

  d = hom(res, A1)
  str = "$(d)"
  str = "$(Oscar.underlying_complex(d))"

  t = tensor_product(res, res, d)
  str = "$(t)"

end
