@testset "Toric schemes" begin
  C = positive_hull([-1 1; 1 1])
  antv = affine_normal_toric_variety(C)

  @testset "A simplicial (and not smooth) affine toric scheme" begin
    @test is_smooth(Oscar.underlying_scheme(antv)) == is_smooth(antv)
    @test dim(Oscar.underlying_scheme(antv)) == dim(antv)
    @test polarize(cone(antv)) == weight_cone(antv)
    @test ngens(ambient_coordinate_ring(antv)) == nrows(hilbert_basis(antv)) == 3
  end

  X = hirzebruch_surface(NormalToricVariety, 3)

  @testset "Toric Scheme of Hirzebruch surface F3" begin
    @test is_smooth(Oscar.underlying_scheme(X)) == is_smooth(X)
    @test dim(Oscar.underlying_scheme(X)) == dim(X)
  end

  IP1 = projective_space(NormalToricVariety, 1)
  set_coordinate_names(IP1, [:x, :y])
  Y = IP1 * IP1

  @testset "Product of projective spaces" begin
    @test is_smooth(Oscar.underlying_scheme(Y)) == is_smooth(Y)
    @test length(values(gluings(default_covering(Y)))) == 16
  end

  IP2 = projective_space(NormalToricVariety, 2)
  set_coordinate_names(IP2, [:x, :y, :z])
  X, iso = Oscar.forget_toric_structure(IP2)

  @testset "Forget toric structure" begin
    @test X isa CoveredScheme
    @test codomain(iso) === IP2
    @test is_smooth(X) == is_smooth(IP2)
    @test length(keys(gluings(default_covering(X)))) == 9
    @test domain(iso) === X
    @test codomain(iso) === IP2
    @test domain(inverse(iso)) === IP2
    @test codomain(inverse(iso)) === X
  end

  IP3 = projective_space(NormalToricVariety, 3)

  @testset "Toric ideal sheaves" begin
    S = cox_ring(IP3)
    I = ideal(S, gens(S)[2:3])
    II = IdealSheaf(IP3, I)
  end

  S = cox_ring(IP2)
  (x, y, z) = gens(S)
  I = IdealSheaf(IP2, ideal(S, [z * (x - y)]))
  J = IdealSheaf(IP2, ideal(S, [x, y]))
  bl = blow_up(IP2, [1, 1])
  pb_I = pullback(bl, I)
  pb_J = pullback(bl, J)

  @testset "Toric blowup morphism as morphism of covered schemes" begin
    @test scheme(I) === IP2
    @test length(Oscar.maximal_associated_points(I)) == 2
    @test length(Oscar.maximal_associated_points(pb_I)) == 3
    @test dim(J) == 0
    @test dim(pb_J) == 1
  end

  F2 = hirzebruch_surface(NormalToricVariety, 2)
  f = toric_morphism(F2, matrix(ZZ, [[1], [0]]), IP1; check=true)
  S = cox_ring(IP1)
  (x, y) = gens(S)
  K = IdealSheaf(IP1, ideal(S, [x]))
  pb_K = pullback(f, K)
  # The support of pb_K should be nothing but P1, it is the fiber
  # of the P1-fibration that defines the Hirebruch surface over the
  # point defined by K.
  # TODO: Add a test that verifies that indeed pb_K is P1.
  @testset "Hirzebruch surface as P1 fibration over P1" begin
    @test length(Oscar.maximal_associated_points(K)) == 1
    @test dim(K) == 0
    @test length(Oscar.maximal_associated_points(pb_K)) == 1
    @test dim(pb_K) == 1
    @test is_smooth(subscheme(pb_K)) == true
  end

  S = cox_ring(IP2)
  x, y, z = gens(S)
  I = ideal(S, [x - y]) * ideal(S, [z])
  II = IdealSheaf(IP2, I)

  @testset "Blowups that leave the toric setting" begin
    @test is_surjective(lattice_homomorphism(Oscar.underlying_morphism(bl)))
    @test is_injective(lattice_homomorphism(Oscar.underlying_morphism(bl)))
    @test length(Oscar.maximal_associated_points(pullback(bl, II))) == 3
    @test length(Oscar.maximal_associated_points(strict_transform(bl, II))) == 2
  end
end

@testset "Lazy gluings" begin
  f = polyhedral_fan(
    incidence_matrix([[1, 2, 3], [1, 4, 5, 6]]),
    [0 0 1; 1 0 1; 0 1 0; -1 0 1; -1 -1 1; 0 -1 1],
  )
  number_of_maximal_cones(f)
  ntv = normal_toric_variety(f)
  X = Oscar.underlying_scheme(ntv)
  for g in values(gluings(default_covering(X)))
    @test !(g isa Oscar.LazyGluing) || !isdefined(g, :G)
  end

  for g in values(gluings(default_covering(X)))
    g isa Oscar.LazyGluing || continue
    g_sub = Oscar.underlying_gluing(g)
    @test g_sub isa Oscar.SimpleGluing
    @test g_sub === g.G
  end
end

@testset "toric divisors to weil divisors" begin
  IP = weighted_projective_space(NormalToricVariety, [3, 2, 5])
  w = canonical_divisor(IP)
  D = Oscar.underlying_divisor(w; check=true)
  D2 = Oscar.underlying_divisor(w; algorithm=:via_polymake, check=true)
  D3 = Oscar.underlying_divisor(w; algorithm=:via_oscar, check=true)
  @test D == D2 == D3

  @test w == forget_toric_structure(w)
  @test w + D == 2 * Oscar.underlying_divisor(w)

  prim = Oscar._torusinvariant_weil_divisors(IP; check=true)
  # Delete the cache manually
  delete!(IP.__attrs, :_torusinvariant_weil_divisors)
  prim2 = Oscar._torusinvariant_weil_divisors(IP; check=true, algorithm=:via_oscar)
  @test prim == prim2
  @test prim[1] !== prim2[1]

  D = WeilDivisor(Oscar._ideal_sheaf_via_polymake(IP, ZZ.([2, 3, 7])))
  D2 = 2 * prim[1] + 3 * prim[2] + 7 * prim[3]

  @test is_effective(D)
  @test is_effective(D2)

  #@test D == D2 # Test takes too long
  IP = projective_space(NormalToricVariety, 1)
  w = canonical_divisor(IP)
  K0 = Oscar.underlying_divisor(w; algorithm=:via_polymake)
  delete!(w.__attrs, :underlying_divisor)
  K1 = Oscar.underlying_divisor(w; algorithm=:direct)
  delete!(w.__attrs, :underlying_divisor)
  K2 = Oscar.underlying_divisor(w; algorithm=:via_oscar)
  delete!(w.__attrs, :underlying_divisor)

  @test K1 == K2
  @test K1 == K0
  @test K2 == K0

  @test w == Oscar.underlying_divisor(w)
end

@testset "general ideal sheaves from ideals in cox ring" begin
  c = cone([1 1; -1 1])
  ntv = normal_toric_variety(c)
  x1, x2 = gens(cox_ring(ntv))
  my_ideal = ideal([x1 * x2])
  is_complete(ntv)
  has_torusfactor(ntv)
  ideal_sheaf(ntv, my_ideal)

  p231 = weighted_projective_space(NormalToricVariety, [2, 3, 1])
  x1, x2, x3 = gens(cox_ring(p231))
  my_ideal = ideal([x1 * x2])
  II = ideal_sheaf(p231, my_ideal)
  JJ = Oscar.ToricIdealSheafFromCoxRingIdeal(p231, my_ideal)
  @test II == JJ

  # The coordinates below correspond to the ideal `ideal([x1, x2])`
  coords = [1, 1, 0]

  pr = blow_up_along_minimal_supercone_coordinates(p231, [1, 1, 0])
  pullback(pr, II)
  @test is_subset(total_transform(pr, II), strict_transform(pr, II))

  X = p231 * p231
  a, b, c, x, y, z = gens(cox_ring(X))
  I = ideal(cox_ring(X), [a^3 + b^2 + c^2 * a^2, a * x, a * y^2 - 25 * c^2 * z^6])
  II = IdealSheaf(X, I)
  JJ = Oscar.ToricIdealSheafFromCoxRingIdeal(X, I)
  @test II == JJ
end

# Test chart restrictions against independent local equations.
@testset "Cox-ring ideal sheaves on toric varieties" begin
  X = projective_space(NormalToricVariety, 2)
  S = cox_ring(X)
  x, y, z = gens(S)
  I = ideal(S, x^2 - y * z)
  II = ideal_sheaf(X, I)
  @test II isa Oscar.ToricIdealSheafFromCoxRingIdeal
  @test Oscar.ideal_in_cox_ring(II) === I

  charts = affine_charts(X)
  @test length(charts) == 3
  for (k, U) in enumerate(charts)
    u, v = gens(OO(U))
    f = k == 2 ? u * v - one(OO(U)) : u^2 - v
    @test II(U) == ideal(OO(U), f)
  end

  zero_sheaf = ideal_sheaf(X, ideal(S, zero(S)))
  unit_sheaf = ideal_sheaf(X, ideal(S, one(S)))
  for U in charts
    @test iszero(zero_sheaf(U))
    @test isone(unit_sheaf(U))
  end

  # The affine quadratic cone is A^2 / {±1}. Its closed torus-fixed point
  # has maximal ideal generated by the three invariant quadratic monomials.
  Y = normal_toric_variety(Oscar.positive_hull([1 1; -1 1]))
  @test !is_smooth(Y)
  T = cox_ring(Y)
  a, b = gens(T)
  point_sheaf = ideal_sheaf(Y, ideal(T, [a, b]))
  V = only(affine_charts(Y))
  @test point_sheaf(V) == ideal(OO(V), gens(OO(V)))

  # AffineNormalToricVariety is a separate, currently unsupported type.
  Y_affine = affine_normal_toric_variety(Oscar.positive_hull([1 1; -1 1]))
  @test_throws MethodError ideal_sheaf(
    Y_affine, ideal(cox_ring(Y_affine), first(gens(cox_ring(Y_affine))))
  )

  # A full-dimensional cone is required for the chart semigroup to be pointed.
  torus_factor = normal_toric_variety(Oscar.positive_hull([1 0]))
  @test !all(is_fulldimensional, maximal_cones(torus_factor))
  @test_throws ArgumentError ideal_sheaf(
    torus_factor, ideal(cox_ring(torus_factor), zero(cox_ring(torus_factor)))
  )

  # Reject ideals from a different ambient Cox ring before lazy evaluation.
  Z = projective_space(NormalToricVariety, 1)
  @test_throws ArgumentError ideal_sheaf(X, ideal(cox_ring(Z), first(gens(cox_ring(Z)))))
end

@testset "Ideal sheaves of toric divisors" begin
  X = projective_space(NormalToricVariety, 2)
  S = cox_ring(X)
  x, y, z = gens(S)
  charts = affine_charts(X)

  # Prime and non-prime divisors are handled uniformly by Cox monomials.
  D = toric_divisor(X, [1, 0, 0])
  J = ideal_sheaf(D)
  @test J isa Oscar.ToricIdealSheafFromCoxRingIdeal
  @test Oscar.ideal_in_cox_ring(J) == ideal(S, x)
  @test sum(isone(J(U)) for U in charts) == 1

  E = toric_divisor(X, [1, 2, 0])
  K = ideal_sheaf(E)
  @test Oscar.ideal_in_cox_ring(K) == ideal(S, x * y^2)
  target = ideal_sheaf(X, ideal(S, x * y^2))
  for U in charts
    @test K(U) == target(U)
  end

  # The zero divisor corresponds to the unit ideal sheaf.
  Z = trivial_divisor(X)
  @test all(isone(ideal_sheaf(Z)(U)) for U in charts)

  @test_throws ArgumentError ideal_sheaf(toric_divisor(X, [-1, 0, 0]))

  # On the quadratic cone, twice a non-Cartier prime divisor is Cartier.
  Y = normal_toric_variety(positive_hull([1 1; -1 1]))
  P = toric_divisor(Y, [1, 0])
  D = toric_divisor(Y, [2, 0])
  @test !is_cartier(P)
  @test is_cartier(D)
  @test_throws ArgumentError ideal_sheaf(P)
  U = only(affine_charts(Y))
  u, v, w = gens(OO(U))
  @test u * v == w^2
  @test ideal_sheaf(D)(U) == ideal(OO(U), u)
end

@testset "Ideal sheaves from toric orbit closures" begin
  X = projective_space(NormalToricVariety, 2)
  S = cox_ring(X)
  Sigma = polyhedral_fan(X)

  tau = first(cones(Sigma, 1))
  @test tau in Sigma
  k = findfirst(==(first(rays(tau))), rays(Sigma))
  I = ideal(S, gens(S)[k])
  J = ideal_sheaf(X, tau)
  @test J isa Oscar.ToricIdealSheafFromCoxRingIdeal
  @test Oscar.ideal_in_cox_ring(J) == I
  expected = ideal_sheaf(X, I)
  for U in affine_charts(X)
    @test J(U) == expected(U)
  end

  sigma = first(maximal_cones(Sigma))
  indices = [findfirst(==(ray), rays(Sigma)) for ray in rays(sigma)]
  K = ideal_sheaf(X, sigma)
  @test Oscar.ideal_in_cox_ring(K) == ideal(S, [gens(S)[i] for i in indices])

  zero_cone = first(cones(Sigma, 0))
  @test zero_cone in Sigma
  @test all(iszero(ideal_sheaf(X, zero_cone)(U)) for U in affine_charts(X))

  invalid_cone = positive_hull([1 2])
  @test !(invalid_cone in Sigma)
  @test_throws ArgumentError ideal_sheaf(X, invalid_cone)

  # The legacy constructor retains its concrete IdealSheaf return type.
  L = IdealSheaf(X, I)
  @test L isa IdealSheaf
  for U in affine_charts(X)
    @test L(U) == J(U)
  end

  # A full-dimensional nonsimplicial cone defines a closed torus-fixed point.
  Y = normal_toric_variety(positive_hull([0 0 1; 1 0 1; 1 1 1; 0 1 1]))
  @test !is_simplicial(Y)
  U = only(affine_charts(Y))
  point_sheaf = ideal_sheaf(Y, first(maximal_cones(Y)))
  @test point_sheaf(U) == ideal(OO(U), gens(OO(U)))
end
