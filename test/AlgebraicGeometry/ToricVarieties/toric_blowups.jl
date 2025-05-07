@testset "Toric blowups" begin
  
  P2 = projective_space(NormalToricVariety, 2)
  BP2 = domain(blow_up(P2, 1; coordinate_name = "e"))
  
  @testset "Basic properties of BP2" begin
    @test is_normal(BP2)
    @test !is_affine(BP2)
    @test is_projective(BP2)
    @test !is_projective_space(BP2) 
    @test is_smooth(BP2)
    @test is_complete(BP2)
    @test is_orbifold(BP2)
    @test is_simplicial(BP2)
    @test !has_torusfactor(BP2)
  end
  
  @testset "Basic attributes of BP2" begin
    @test betti_number(BP2, 0) == 1
    @test betti_number(BP2, 1) == 0
    @test betti_number(BP2, 2) == 2
    @test betti_number(BP2, 3) == 0
    @test betti_number(BP2, 4) == 1
    @test euler_characteristic(BP2) == 4
    @test torsion_free_rank(picard_group(BP2)) == 2
  end

  amb = load(joinpath(Oscar.oscardir, "test/AlgebraicGeometry/ToricVarieties", "pr3006.ntv"))
  (x1, x2, x3, x4, x, y, z) = gens(cox_ring(amb));

  # The coordinates below correspond to the ideal `ideal([x, y, x1]`
  coords = [1, 0, 0, 0, 1, 1, 0]

  bd1 = blow_up_along_minimal_supercone_coordinates(
    amb, coords; coordinate_name = "e1"
  )
  amb1 = domain(bd1)
  
  @testset "Trigger issue of PR3006" begin
    @test is_complete(amb1)
    @test is_simplicial(amb1)
  end
  
  P2 = projective_space(NormalToricVariety, 2)
  bl = blow_up(P2, [1, 1])
  
  @testset "Basic tests for simple toric blowup" begin
    @test torsion_free_rank(domain(lattice_homomorphism(bl))) == 2
    @test torsion_free_rank(codomain(lattice_homomorphism(bl))) == 2
    @test rank(matrix(morphism_on_torusinvariant_weil_divisor_group(bl))) == 3
    @test matrix(morphism_on_torusinvariant_cartier_divisor_group(bl)) == matrix(morphism_on_torusinvariant_weil_divisor_group(bl))
  end
  
  bl2 = blow_up(P2, [-1, 0])

  @testset "Arithmetics, comparison and composition" begin
    @test lattice_homomorphism(bl + bl) == lattice_homomorphism(2 * bl)
    @test lattice_homomorphism(bl - bl) == lattice_homomorphism(0 * bl)
    @test lattice_homomorphism(bl * toric_identity_morphism(codomain(bl))) == lattice_homomorphism(bl)
    @test bl !== bl2
  end

  S = cox_ring(P2)
  J = IdealSheaf(P2, ideal(S, S[1]))
  pbJ = Oscar.total_transform(bl, J)
  pbJ_str = strict_transform(bl, J)
  E = exceptional_prime_divisor(bl)
  H = toric_divisor(domain(bl), [1, 2, 3, 4]) + E
  
  @testset "Strict, total transform and exceptional divisor for simple toric blowup" begin
    @test issubset(ideal_sheaf(H), ideal_sheaf(E))
    @test issubset(pbJ, pbJ_str)
    @test issubset(pbJ, ideal_sheaf(E))
    @test !issubset(pbJ_str, ideal_sheaf(E))
  end

  @testset "Toric blowups along singular cones" begin
    # 1/2(1, 1) quotient singularity, blowup along the maximal cone
    ray_generators = [[2, -1], [0, 1]]
    max_cones = IncidenceMatrix([[1, 2]])
    X = normal_toric_variety(max_cones, ray_generators)
    f = blow_up(X, 1)
    @test ray_vector(QQFieldElem, [1, 0]) in rays(domain(f))
    @test (
      minimal_supercone_coordinates_of_exceptional_ray(f)
      ==
      QQFieldElem[1//2, 1//2]
    )
    
    # Now blowing up along an existing ray
    g = blow_up(X, 2)
    @test n_rays(domain(g)) == 2
    @test n_cones(domain(g)) == 4

    # Quadratic cone, blowup along maximal cone
    ray_generators = [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    max_cones = IncidenceMatrix([[1, 2, 3, 4]])
    PF = polyhedral_fan(max_cones, ray_generators)
    X = normal_toric_variety(PF)
    f = blow_up(X, 1)
    @test ray_vector(QQFieldElem, [1, 1, 2]) in rays(domain(f))
    @test (
      minimal_supercone_coordinates_of_exceptional_ray(f)
      ==
      QQFieldElem[1//2, 1//2, 1//2, 1//2]
    )
    
    # Now blowing up along an existing ray
    g = blow_up(X, 6)
    @test n_rays(domain(g)) == 4
    @test n_cones(domain(g)) == 12
  end
end
