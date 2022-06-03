@testset "Secondary invariants (for matrix groups)" begin
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

    # The secondary invariants have to be a module basis. Let's test this for
    # the degree 9 homogeneous component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, [ f.f for f in primary_invariants(RG) ])
    b1, _ = Oscar.generators_for_given_degree!(C, [ f.f for f in secondary_invariants(RG) ], 9, false)
    b2 = [ f.f for f in basis(RG, 9) ]
    B = Oscar.BasisOfPolynomials(R, b1)
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  s_invars = secondary_invariants(RGm)
  actionN3 = Oscar.right_action(polynomial_ring(RGm), N3)
  actionN4 = Oscar.right_action(polynomial_ring(RGm), N4)
  for f in s_invars
    @test actionN3(f) == f
    @test actionN4(f) == f
  end
  @test length(s_invars) == 6
  # The secondary invariants have to be a module basis. Let's test this for
  # the degree 6 homogeneous component.
  R = polynomial_ring(RGm).R
  C = Oscar.PowerProductCache(R, [ f.f for f in primary_invariants(RGm) ])
  b1, _ = Oscar.generators_for_given_degree!(C, [ f.f for f in secondary_invariants(RGm) ], 6, false)
  b2 = [ f.f for f in basis(RGm, 6) ]
  B = Oscar.BasisOfPolynomials(R, b1)
  for f in b2
    @test !Oscar.add_to_basis!(B, f)
  end

  for RG in [ RG0, RGp ]
    is_invars = irreducible_secondary_invariants(RG)
    @test issubset(is_invars, secondary_invariants(RG))
    @test length(is_invars) == 1

    # The irreducible secondary invariants have to be algebra generators. Let's
    # test this for the degree 9 homogeneous component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, append!([ f.f for f in primary_invariants(RG) ], [ f.f for f in is_invars ]))
    b1 = Oscar.all_power_products_of_degree!(C, 9, false)
    b2 = [ f.f for f in basis(RG, 9) ]
    B = Oscar.BasisOfPolynomials(R, b1)
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  is_invars = irreducible_secondary_invariants(RGm)
  @test issubset(is_invars, secondary_invariants(RGm))
  @test length(is_invars) == 2
  # The irreducible secondary invariants have to be algebra generators. Let's
  # test this for the degree 9 homogeneous component.
  R = polynomial_ring(RGm).R
  C = Oscar.PowerProductCache(R, append!([ f.f for f in primary_invariants(RGm) ], [ f.f for f in is_invars ]))
  b1 = Oscar.all_power_products_of_degree!(C, 6, false)
  b2 = [ f.f for f in basis(RGm, 6) ]
  B = Oscar.BasisOfPolynomials(R, b1)
  for f in b2
    @test !Oscar.add_to_basis!(B, f)
  end

  for RG in [ RG0, RGp, RGm ]
    is_invars = irreducible_secondary_invariants(RG)
    for i = 1:length(RG.secondary.invars)
      f = RG.secondary.invars[i]
      @test in(f, is_invars) == RG.secondary.is_irreducible[i]
      t = one(polynomial_ring(RG))
      for j = 1:length(is_invars)
        t *= is_invars[j]^RG.secondary.sec_in_irred[i][j]
      end
      @test f == t
    end
  end

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

@testset "Secondary invariants (for permutation groups)" begin
  G = sylow_subgroup(symmetric_group(6), 2)[1]
  RG0 = invariant_ring(G)     # char. 0

  F7 = GF(7)
  RGp = invariant_ring(F7, G) # char. p, non-modular

  F9 = GF(3, 2)
  RGm = invariant_ring(F9, G) # char. p, modular

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

    # The secondary invariants have to be a module basis. Let's test this for
    # the degree 9 homogeneous component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, [ f.f for f in primary_invariants(RG) ])
    b1, _ = Oscar.generators_for_given_degree!(C, [ f.f for f in secondary_invariants(RG) ], 9, false)
    b2 = [ f.f for f in basis(RG, 9) ]
    B = Oscar.BasisOfPolynomials(R, b1)
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  s_invars = secondary_invariants(RGm)
  actions = [Oscar.right_action(polynomial_ring(RGm), x) for x in gens(G)]
  for f in s_invars
    @test all(act -> act(f) == f, actions)
  end
  @test length(s_invars) == 2
  # The secondary invariants have to be a module basis. Let's test this for
  # the degree 9 homogeneous component.
  R = polynomial_ring(RGm).R
  C = Oscar.PowerProductCache(R, [ f.f for f in primary_invariants(RGm) ])
  b1, _ = Oscar.generators_for_given_degree!(C, [ f.f for f in secondary_invariants(RGm) ], 9, false)
  b2 = [ f.f for f in basis(RGm, 9) ]
  B = Oscar.BasisOfPolynomials(R, b1)
  for f in b2
    @test !Oscar.add_to_basis!(B, f)
  end

  for RG in [ RG0, RGp ]
    is_invars = irreducible_secondary_invariants(RG)
    @test issubset(is_invars, secondary_invariants(RG))
    @test length(is_invars) == 1

    # The irreducible secondary invariants have to be algebra generators. Let's
    # test this for the degree 9 homogeneous component.
    R = polynomial_ring(RG).R
    C = Oscar.PowerProductCache(R, append!([ f.f for f in primary_invariants(RG) ], [ f.f for f in is_invars ]))
    b1 = Oscar.all_power_products_of_degree!(C, 9, false)
    b2 = [ f.f for f in basis(RG, 9) ]
    B = Oscar.BasisOfPolynomials(R, b1)
    for f in b2
      @test !Oscar.add_to_basis!(B, f)
    end
  end

  is_invars = irreducible_secondary_invariants(RGm)
  @test issubset(is_invars, secondary_invariants(RGm))
  @test length(is_invars) == 1
  # The irreducible secondary invariants have to be algebra generators. Let's
  # test this for the degree 9 homogeneous component.
  R = polynomial_ring(RGm).R
  C = Oscar.PowerProductCache(R, append!([ f.f for f in primary_invariants(RGm) ], [ f.f for f in is_invars ]))
  b1 = Oscar.all_power_products_of_degree!(C, 6, false)
  b2 = [ f.f for f in basis(RGm, 6) ]
  B = Oscar.BasisOfPolynomials(R, b1)
  for f in b2
    @test !Oscar.add_to_basis!(B, f)
  end

  for RG in [ RG0, RGp, RGm ]
    is_invars = irreducible_secondary_invariants(RG)
    for i = 1:length(RG.secondary.invars)
      f = RG.secondary.invars[i]
      @test in(f, is_invars) == RG.secondary.is_irreducible[i]
      t = one(polynomial_ring(RG))
      for j = 1:length(is_invars)
        t *= is_invars[j]^RG.secondary.sec_in_irred[i][j]
      end
      @test f == t
    end
  end

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
  ggens, exps = Oscar.generators_for_given_degree!(C, [ x[3] ], 3, true)
  @test Set(ggens) == Set([ x[1]^2*x[3], x[1]*x[2]*x[3], x[2]^2*x[3] ])
  for m in ggens
    @test haskey(exps, m)
    @test set_exponent_vector!(one(R), 1, exps[m]) == m
  end
end
