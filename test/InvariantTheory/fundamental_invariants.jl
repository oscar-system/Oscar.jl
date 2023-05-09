@testset "Fundamental invariants (for matrix groups)" begin
  # Char 0
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 0, 0, 1, 1, 0, 0 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])
  RG0 = invariant_ring(M1, M2)

  # Call it once without specifying `algo`
  @test length(fundamental_invariants(RG0)) == 4

  for algo in [ :king, :primary_and_secondary ]
    RG = invariant_ring(M1, M2) # redefine to avoid caching

    invars = fundamental_invariants(RG, algo)
    @test length(invars) == 4
    @test [ total_degree(f.f) for f in invars ] == [ 3, 3, 6, 9 ]
    for f in invars
      @test reynolds_operator(RG, f) == f
    end

    # Let's test whether invars are algebra generators for the degree 9 homogeneous
    # component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, [ f.f for f in invars ])
    B = Oscar.BasisOfPolynomials(R, Oscar.all_power_products_of_degree!(C, 9, false))
    b2 = [ f.f for f in basis(RG, 9) ]
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  # Char p, non-modular
  F = GF(3)
  N1 = matrix(F, 3, 3, [ 0, 1, 0, 2, 0, 0, 0, 0, 2 ])
  N2 = matrix(F, 3, 3, [ 2, 0, 0, 0, 2, 0, 0, 0, 2 ])

  for algo in [ :king, :primary_and_secondary ]
    RG = invariant_ring(N1, N2) # redefine to avoid caching

    invars = fundamental_invariants(RG, algo)
    @test length(invars) == 4
    @test [ total_degree(f.f) for f in invars ] == [ 2, 2, 4, 4 ]
    for f in invars
      @test reynolds_operator(RG, f) == f
    end

    # Let's test whether invars are algebra generators for the degree 6 homogeneous
    # component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, [ f.f for f in invars ])
    B = Oscar.BasisOfPolynomials(R, Oscar.all_power_products_of_degree!(C, 6, false))
    b2 = [ f.f for f in basis(RG, 6) ]
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  # Char p, modular
  F9, b = FiniteField(3, 2, "b")
  N3 = matrix(F9, [ 1 0 0 0; b + 1 1 0 0; -1 0 1 0; b 0 -1 1 ])
  N4 = matrix(F9, [ 1 0 0 0; 1 1 0 0; 1 0 1 0; b -b b 1 ])
  RGm = invariant_ring(N3, N4)

  @test_throws AssertionError fundamental_invariants(RGm, :king)

  invars = fundamental_invariants(RGm)
  @test length(invars) == 6
  @test [ total_degree(f.f) for f in invars ] == [ 1, 2, 3, 3, 4, 9 ]

  actionN3 = Oscar.right_action(polynomial_ring(RGm), N3)
  actionN4 = Oscar.right_action(polynomial_ring(RGm), N4)
  for f in invars
    @test actionN3(f) == f
    @test actionN4(f) == f
  end

  # Let's test whether invars are algebra generators for the degree 6 homogeneous
  # component.
  R = polynomial_ring(RGm).R
  C = Oscar.PowerProductCache(R, [ f.f for f in invars ])
  B = Oscar.BasisOfPolynomials(R, Oscar.all_power_products_of_degree!(C, 6, false))
  b2 = [ f.f for f in basis(RGm, 6) ]
  for f in b2
    @test !Oscar.add_to_basis!(B, f)
  end

  # Test some special cases
  # Cyclic group in King's algorithm
  M = matrix(QQ, [ 0 1 ; 1 0 ])
  RG = invariant_ring(M)

  invars = fundamental_invariants(RG, :king)
  @test length(invars) == 2
  @test [ total_degree(f.f) for f in invars ] == [ 1, 2 ]
  for f in invars
    @test reynolds_operator(RG, f) == f
  end

  # Specify degree bound
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])

  invars1 = fundamental_invariants(invariant_ring(M1, M2))
  invars2 = fundamental_invariants(invariant_ring(M1, M2), beta = 6)

  @test [ f(gens(parent(invars1[1]))...) for f in invars2 ] == invars1

  # No irreducible secondary invariants
  M = matrix(QQ, [ 0 1 ; 1 0 ])
  RG = invariant_ring(M)

  invars = fundamental_invariants(RG, :primary_and_secondary)
  @test length(invars) == 2
  @test [ total_degree(f.f) for f in invars ] == [ 1, 2 ]
  for f in invars
    @test reynolds_operator(RG, f) == f
  end

  # `:primary_and_secondary` actually removes invariants
  K, a = CyclotomicField(3, "a")
  M1 = matrix(K, 3, 3, [ 0, 1, 0, 1, 0, 0, 0, 0, 1 ])
  M2 = matrix(K, 3, 3, [ 1, 0, 0, 0, a, 0, 0, 0, -a - 1 ])

  RG = invariant_ring(M1, M2)
  primary_invariants(RG, primary_degrees = [ 3, 6, 6 ])

  invars = fundamental_invariants(RG, :primary_and_secondary)
  @test [ total_degree(f.f) for f in invars ] == [ 3, 3, 3, 6 ]

  S = RG.fundamental.S
  RtoS = RG.fundamental.toS
  R = polynomial_ring(RG)
  StoR = hom(S, R, invars)
  for (f, g) in RtoS
    @test StoR(g) == f
  end
end

@testset "Fundamental invariants (for permutation groups)" begin
  G = sylow_subgroup(symmetric_group(6), 2)[1]

  # Char 0
  for algo in [ :king, :primary_and_secondary ]
    RG = invariant_ring(G) # redefine to avoid caching

    invars = fundamental_invariants(RG, algo)
    @test length(invars) == 7
    @test [ total_degree(f.f) for f in invars ] == [ 1, 1, 2, 2, 2, 3, 4 ]
    for f in invars
      @test reynolds_operator(RG, f) == f
    end

    # Let's test whether invars are algebra generators for the degree 9 homogeneous
    # component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, [ f.f for f in invars ])
    B = Oscar.BasisOfPolynomials(R, Oscar.all_power_products_of_degree!(C, 9, false))
    b2 = [ f.f for f in basis(RG, 9) ]
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  # Char p, non-modular
  for algo in [ :king, :primary_and_secondary ]
    RG = invariant_ring(GF(7), G) # redefine to avoid caching

    invars = fundamental_invariants(RG, algo)
    @test length(invars) == 7
    @test [ total_degree(f.f) for f in invars ] == [ 1, 1, 2, 2, 2, 3, 4 ]
    for f in invars
      @test reynolds_operator(RG, f) == f
    end

    # Let's test whether invars are algebra generators for the degree 9 homogeneous
    # component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, [ f.f for f in invars ])
    B = Oscar.BasisOfPolynomials(R, Oscar.all_power_products_of_degree!(C, 9, false))
    b2 = [ f.f for f in basis(RG, 9) ]
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  # Char p, modular
  RG = invariant_ring(GF(2), G)
  @test_throws AssertionError fundamental_invariants(RG, :king)

  invars = fundamental_invariants(RG)
  @test length(invars) == 7
  @test [ total_degree(f.f) for f in invars ] == [ 1, 1, 2, 2, 2, 3, 4 ]

  actions = [ Oscar.right_action(polynomial_ring(RG), x) for x in gens(G) ]
  for f in invars
    @test all(act -> act(f) == f, actions)
  end

  # Let's test whether invars are algebra generators for the degree 9 homogeneous
  # component.
  R = polynomial_ring(RG).R
  C = Oscar.PowerProductCache(R, [ f.f for f in invars ])
  B = Oscar.BasisOfPolynomials(R, Oscar.all_power_products_of_degree!(C, 9, false))
  b2 = [ f.f for f in basis(RG, 9) ]
  for f in b2
    @test !Oscar.add_to_basis!(B, f)
  end
end
