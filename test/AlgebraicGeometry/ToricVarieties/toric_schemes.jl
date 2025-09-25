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
