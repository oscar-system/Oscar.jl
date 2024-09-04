@testset "Toric blowdowns" begin
  
  P2 = projective_space(NormalToricVariety, 2)
  BP2 = domain(blow_up(P2, 1; coordinate_name = "e"))
  
  @testset "Basic properties of BP2" begin
    @test is_normal(BP2) == true
    @test is_affine(BP2) == false
    @test is_projective(BP2) == true
    @test is_projective_space(BP2) == false
    @test is_smooth(BP2) == true
    @test is_complete(BP2) == true
    @test is_orbifold(BP2) == true
    @test is_simplicial(BP2) == true
    @test has_torusfactor(BP2) == false
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
  bd1 = blow_up(amb, ideal([x, y, x1]); coordinate_name = "e1")
  amb1 = domain(bd1)
  
  @testset "Trigger issue of PR3006" begin
    @test is_complete(amb1)
    @test is_simplicial(amb1)
  end
  
  P2 = projective_space(NormalToricVariety, 2)
  bl = blow_up(P2, [1, 1])
  
  @testset "Basic tests for simple toric blowdown" begin
    @test torsion_free_rank(domain(grid_morphism(bl))) == 2
    @test torsion_free_rank(codomain(grid_morphism(bl))) == 2
    @test matrix(morphism_on_torusinvariant_weil_divisor_group(bl)) == matrix(ZZ, [1 0 0; 0 1 0; 0 0 1; 1 1 0])
    @test matrix(morphism_on_torusinvariant_cartier_divisor_group(bl)) == matrix(morphism_on_torusinvariant_weil_divisor_group(bl))
  end
  
  bl2 = blow_up(P2, [-1, 0])

  @testset "Arithmetics, comparison and composition" begin
    @test grid_morphism(bl + bl) == grid_morphism(2 * bl)
    @test grid_morphism(bl - bl) == grid_morphism(0 * bl)
    @test grid_morphism(bl * toric_identity_morphism(codomain(bl))) == grid_morphism(bl)
    @test bl !== bl2
  end

  S = cox_ring(P2)
  J = IdealSheaf(P2, ideal(S, S[1]))
  pbJ = Oscar.total_transform(bl, J)
  pbJ_str = strict_transform(bl, J)
  E = exceptional_divisor(bl)
  H = toric_divisor(domain(bl), [1, 2, 3, 4]) + E
  
  @testset "Strict, total transform and exceptional divisor for simple toric blowdown" begin
    @test issubset(ideal_sheaf(H), ideal_sheaf(E))
    @test issubset(pbJ, pbJ_str)
    @test issubset(pbJ, ideal_sheaf(E))
    @test !issubset(pbJ_str, ideal_sheaf(E))
  end

  I = ideal(S, [S[1], S[2]])
  bl2 = blow_up(P2, I)
  II = IdealSheaf(P2, I)
  
  @testset "Properties of toric blowdown defined by ideal" begin
    @test II == center(bl)
    @test bl2 isa Oscar.ToricBlowdownMorphism
    @test center(bl2) == II
  end
end
