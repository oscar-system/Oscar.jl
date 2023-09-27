@testset "Blowup of toric varieties (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    P2 = projective_space(NormalToricVariety, 2; set_attributes)
    BP2 = domain(blow_up(P2, 1; coordinate_name = "e", set_attributes = set_attributes))
    
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
        @test rank(picard_group(BP2)) == 2
    end
end

@testset "Extended functionality for ToricBlowdownMorphism" begin
    IP2 = projective_space(NormalToricVariety, 2)
    bl = blow_up(IP2, [1, 1])

    grid_morphism(bl)
    domain(bl)
    codomain(bl)
    morphism_on_torusinvariant_weil_divisor_group(bl)
    morphism_on_torusinvariant_cartier_divisor_group(bl)
    covering_morphism(bl)
    exceptional_divisor(bl)
    
    S = cox_ring(IP2)
    I = ideal(S, [S[1], S[2]])
    II = IdealSheaf(IP2, I)
    @test II == center(bl)
    
    bl2 = blow_up(IP2, I)
    @test bl2 isa Oscar.ToricBlowdownMorphism
    @test center(bl2) == II

    E = exceptional_divisor(bl)
    X = domain(bl)
    ideal_sheaf(E)
    D = toric_divisor(X, [1, 2, 3, 4])
    H = D + E
    @test issubset(ideal_sheaf(H), ideal_sheaf(E))
    J = IdealSheaf(IP2, ideal(S, S[1]))
    pbJ = Oscar.total_transform(bl, J)
    pbJ_str = strict_transform(bl, J)
    @test issubset(pbJ, pbJ_str)
    @test issubset(pbJ, ideal_sheaf(E))
    @test !issubset(pbJ_str, ideal_sheaf(E))
end
