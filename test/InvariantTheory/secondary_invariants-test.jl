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
end
