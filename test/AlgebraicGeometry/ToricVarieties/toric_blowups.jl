@testset "Blowup of toric varieties" begin
  
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

  @testset "Trigger issue of PR3006" begin
    amb = load(joinpath(Oscar.oscardir, "test/AlgebraicGeometry/ToricVarieties", "pr3006.ntv"))
    (x1, x2, x3, x4, x, y, z) = gens(cox_ring(amb));
    bd1 = blow_up(amb, ideal([x, y, x1]); coordinate_name = "e1")
    amb1 = domain(bd1)
    @test is_complete(amb1)
    @test is_simplicial(amb1)
  end
end
