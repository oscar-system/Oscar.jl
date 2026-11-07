using Test
using Oscar

@testset "native primary decomposition" begin
  # Test the native implementation of the algorithms of `primdec.lib`
  # (GTZ, SY, radical and equidimensional decompositions) against Singular.
  same_ideals(A, B) = length(A) == length(B) && all(a -> any(b -> a == b, B), A)
  for K in (QQ, GF(32003))
    R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
    examples = [
      ideal(R, [x^2, x*y]),
      ideal(R, [x + y + z, x*y + y*z + z*x, x*y*z - 1]),
      ideal(R, [x*(x-1)*(x-2), y - x^2, z - x*y]),
      ideal(R, [x^2*y, x*y^2, x*z^2, z^3*y]),
      ideal(R, [x^2 - y*z, y^3 - z^2*x, x*y*z]),
      ideal(R, [(x^2+y^2+z^2-1)*(x-y), (x^2+y^2+z^2-1)*z]),
      ideal(R, [(x^2-2)*(y^2+1), (x^2-2)*z, (y^2+1)*z^2]),
    ]
    for I in examples
      ref = primary_decomposition(I)
      ref_primes = [P for (_, P) in ref]
      @test same_ideals(minimal_associated_primes_gtz(I), minimal_primes(I))
      for alg in (:gtz, :sy)
        pd = primary_decomposition_native(I; algorithm=alg)
        @test length(pd) == length(ref)
        @test same_ideals([P for (_, P) in pd], ref_primes)
        @test Oscar.test_primary(pd, I)
      end
      @test Oscar.radical_kl(I) == radical(I)
      @test Oscar.radical_equidim(I) == radical(I)
    end
  end

  # equidimensional decompositions
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [x*y, x*z])
  @test length(Oscar.equidim(I)) == 2
  @test Oscar.equidim_max(I) == ideal(R, [x])
  @test Oscar.equi_radical(I) == ideal(R, [x])
  @test intersect(Oscar.prepare_ass(I)...) == radical(I)
  @test Oscar.equidim(ideal(R, [one(R)])) == [ideal(R, [one(R)])]

  # unit ideal, the zero ideal and zero-dimensional ideals
  @test isempty(primary_decomposition_gtz(ideal(R, [one(R)])))
  @test isempty(minimal_associated_primes_gtz(ideal(R, [one(R)])))
  @test length(minimal_associated_primes_gtz(ideal(R, elem_type(R)[]))) == 1
  J = ideal(R, [(x^2-2)*(x^2-3)^2, y^2 - 2 - x, z])
  @test length(primary_decomposition_gtz(J)) == 2
  @test Oscar.test_primary(primary_decomposition_sy(J), J)
  @test_throws ArgumentError primary_decomposition_native(I; algorithm=:foo)
end
