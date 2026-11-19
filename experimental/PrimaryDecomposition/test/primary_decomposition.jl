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

  # no spurious lower-dimensional loci in the equidimensional decomposition
  R4, (a, b, c, d) = polynomial_ring(QQ, [:a, :b, :c, :d])
  I4 = intersect(ideal(R4, [a, b]), ideal(R4, [b, c])^2, ideal(R4, [a - c, b - 1]))
  E4 = Oscar.equidim(I4)
  @test length(E4) == 1
  @test intersect(E4...) == I4

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

  # the nilradical from the trace form agrees with the one from the radical generators
  for K in (QQ, GF(32003))
    R, (x, y) = polynomial_ring(K, [:x, :y])
    I = ideal(R, [x^2 * (x^2 - 2), (y^2 - x)^2])
    Z = Oscar._pd_zero_data(I)
    T, _ = polynomial_ring(K, :T; cached=false)
    Nt = Oscar._pd_nilradical_trace(Z)
    Nr = Oscar._pd_nilradical(Z, Oscar._pd_radical_generators(Z, T))
    @test Nt !== nothing && Nr !== nothing
    @test Nt[1] == Nr[1]
    @test Nt[2][1:Nt[1], :] == Nr[2][1:Nr[1], :]
  end

  # Singular's zero_decomp on a lexicographical basis agrees with the native decomposition
  for K in (QQ, GF(32003))
    R, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
    I = intersect(ideal(R, [x^2 - y*z, y^2 - z*w^2]), ideal(R, [x - w, y + z]), ideal(R, [z^2 + 1, x - y]))
    U = Oscar._pd_independent_set(I)
    E = Oscar._pd_extension(R, U)
    new, _ = Oscar._pd_zero_decomposition_lex(E, I; primary=true)
    old = Oscar._pd_zero_dimensional_decomposition(Oscar._pd_extend_ideal(E, I))
    Pn = [Oscar._pd_contract_lex(E, P) for (_, P) in new]
    Po = [Oscar._pd_contract(E, gens(P)) for (_, P) in old]
    Qn = [Oscar._pd_contract_lex(E, Q) for (Q, _) in new]
    Qo = [Oscar._pd_contract(E, gens(Q)) for (Q, _) in old]
    @test length(Pn) == length(Po)
    @test all(P -> any(==(P), Po), Pn)
    @test all(Q -> any(==(Q), Qo), Qn)
  end

  # contraction from generators found by linear algebra agrees with saturation
  R, (x, y, z, w) = polynomial_ring(GF(32003), [:x, :y, :z, :w])
  I = intersect(ideal(R, [x^2 - y*z, y^2 - z*w^2]), ideal(R, [x - w, y + z]))
  U = Oscar._pd_independent_set(I)
  E = Oscar._pd_extension(R, U)
  Jx = Oscar._pd_extend_ideal(E, I)
  Jr = Oscar._pd_zero_radical(Jx)
  L = ideal(R, [Oscar._pd_from_ext(E, g) for g in gens(Jr)])
  S = Oscar._pd_saturate_leads(L, U)
  new = Oscar._pd_contraction_generators(E, Oscar._pd_zero_data(Jr), maximum(total_degree, gens(L)))
  @test all(g -> g in S, new)
  @test Oscar._pd_saturate_leads(L + ideal(R, new), U) == S
  @test Oscar._pd_contract(E, gens(Jr)) == S

  # the primes of the lex decomposition are radicals of the components, and a
  # decomposition that loses components is detected: the intersection of powers of
  # four linear primes has these four minimal primes, and cyclic-6 has 156 solutions
  R, x = polynomial_ring(GF(32003), 6, :x)
  Ps = [ideal(R, [x[1], x[2]]), ideal(R, [x[3] - x[1], x[4] + x[2]]),
        ideal(R, [x[5] - x[1] - x[3], x[6] - x[2]]), ideal(R, [x[1] + x[2] + x[3], x[4] - x[5], x[6]])]
  I = intersect(Ps[1]^3, Ps[2]^3, Ps[3]^2, Ps[4]^2)
  L = minimal_associated_primes_gtz(I)
  @test length(L) == 4
  @test all(P -> any(==(P), L), Ps)
  R, x = polynomial_ring(GF(32003), 6, :x)
  I = ideal(R, push!([sum(prod(x[mod1(i + j, 6)] for j in 0:k-1) for i in 1:6) for k in 1:5], prod(x) - 1))
  L = minimal_associated_primes_gtz(I)
  @test length(L) == 90
  @test sum(degree(P) for P in L) == 156
end
