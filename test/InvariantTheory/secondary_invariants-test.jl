@testset "Secondary invariants" begin
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])
  RG0 = invariant_ring(M1, M2)

  F3 = GF(3)
  N1 = matrix(F3, 3, 3, [ 0, 1, 0, 2, 0, 0, 0, 0, 2 ])
  N2 = matrix(F3, 3, 3, [ 2, 0, 0, 0, 2, 0, 0, 0, 2 ])
  RGp = invariant_ring(N1, N2) # char p, non-modular

  F9, b = FiniteField(3, 2, "b")
  N3 = matrix(F9, [ 1 0 0 0; b + 1 1 0 0; -1 0 1 0; b 0 -1 1 ])
  N4 = matrix(F9, [ 1 0 0 0; 1 1 0 0; 1 0 1 0; b -b b 1 ])
  RGm = invariant_ring(N3, N4) # char p, modular

  for RG in [ RG0, RGp ]
    s_invars = secondary_invariants(RG)
    m = molien_series(RG)
    S = base_ring(parent(m))
    t = gen(S)
    n = S()
    for f in s_invars
      @test reynolds_operator(RG, f) == f
      n += t^total_degree(f.f)
    end
    d = prod( 1 - t^total_degree(f.f) for f in primary_invariants(RG) )
    @test m == n//d
  end

  s_invars = secondary_invariants(RGm)
  actionN3 = Oscar.right_action(polynomial_ring(RGm), N3)
  actionN4 = Oscar.right_action(polynomial_ring(RGm), N4)
  for f in s_invars
    @test actionN3(f) == f
    @test actionN4(f) == f
  end
  @test length(s_invars) == 6

  is_invars = irreducible_secondary_invariants(RG0)
  @test issubset(is_invars, secondary_invariants(RG0))
  @test length(is_invars) == 1

  is_invars = irreducible_secondary_invariants(RGp)
  @test issubset(is_invars, secondary_invariants(RGp))
  @test length(is_invars) == 1

  is_invars = irreducible_secondary_invariants(RGm)
  @test issubset(is_invars, secondary_invariants(RGm))
  @test length(is_invars) == 2

  R, x = PolynomialRing(QQ, "x" => 1:3)
  C = Oscar.PowerProductCache(R, x)
  @test Set(Oscar.all_power_products_of_degree!(C, 3, false)) == Set(collect(Oscar.all_monomials(R, 3)))
  mons = Oscar.all_power_products_of_degree!(C, 3, true)
  @test Set(mons) == Set(collect(Oscar.all_monomials(R, 3)))
  for m in mons
    @test haskey(C.exponent_vectors, m)
    @test set_exponent_vector!(one(R), 1, C.exponent_vectors[m]) == m
  end

  C = Oscar.PowerProductCache(R, [ x[1], x[2] ])
  gens, exps = Oscar.generators_for_given_degree!(C, [ x[3] ], 3, true)
  @test Set(gens) == Set([ x[1]^2*x[3], x[1]*x[2]*x[3], x[2]^2*x[3] ])
  for m in gens
    @test haskey(exps, m)
    @test set_exponent_vector!(one(R), 1, exps[m]) == m
  end
end
