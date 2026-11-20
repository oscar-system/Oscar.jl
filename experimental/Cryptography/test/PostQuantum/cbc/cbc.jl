using Random
using Test
using Oscar

@testset "cbc goppa codes" begin
  # Test that Patterson's algorithm corrects up to t errors of a Goppa code
  rng = Random.Xoshiro(1)
  for (m, n, t) in ((4, 16, 2), (5, 32, 3), (6, 64, 5))
    C = Oscar._goppa_code(rng, m, n, t)
    @test size(C.H) == (m * t, n)
    @test size(C.G) == (n - m * t, n)
    @test is_zero(C.G * transpose(C.H))
    for w in 0:t
      pos = sort(randperm(rng, n)[1:w])
      e = zeros(Int, n)
      e[pos] .= 1
      s = [is_one(c) ? 1 : 0 for c in vec(collect(C.H * matrix(GF(2), n, 1, e)))]
      @test Oscar._goppa_decode(C, s) == pos
    end
  end
end
