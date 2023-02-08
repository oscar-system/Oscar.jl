@testset "(Weighted) projective space" begin

    P2 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
    P2v2 = projective_space(NormalToricVariety, 2)
    P3 = projective_space(NormalToricVariety, 3)

    @testset "P2" begin
        @test is_normal(P2) == true
        @test is_affine(P2) == false
        @test is_projective(P2) == true
        @test is_smooth(P2) == true
        @test is_complete(P2) == true
        @test has_torusfactor(P2) == false
        @test is_orbifold(P2) == true
        @test is_simplicial(P2) == true
        @test betti_number(P2, 0) == 1
        @test betti_number(P2, 1) == 0
        @test betti_number(P2, 2) == 1
        @test betti_number(P2, 3) == 0
        @test betti_number(P2, 4) == 1
        @test ngens(cox_ring(P2)) == 3
        @test length(stanley_reisner_ideal(P2).gens) == 1
        @test length(irrelevant_ideal(P2).gens) == 3
        @test is_projective_space(P2) == true
    end

    @testset "Weighted projective space" begin
        WPS = weighted_projective_space(NormalToricVariety, [2, 3, 1])
        @test is_smooth(WPS) == false
        @test ngens(cox_ring(WPS)) == 3
        for d in ([1,2,3],[3,1,2],[3,45,68,7,15],[3,1,2,6,5,1,1,8])
            wps = weighted_projective_space(NormalToricVariety, d)
            rm = matrix(ZZ, rays(wps))
            @test wps isa NormalToricVariety
            @test dim(wps) == length(d) - 1
            @test iszero(matrix(ZZ, [d]) * rm)
        end
    end

    @testset "Standard constructors" begin
        @test matrix(map_from_torusinvariant_weil_divisor_group_to_class_group(P2v2)) == matrix(ZZ, [[1], [1], [1]])
        @test transpose(matrix(ZZ,rays(P2v2))) == matrix(map_from_character_lattice_to_torusinvariant_weil_divisor_group(P2v2))
        @test domain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(P2v2)) == character_lattice(P2v2)
        @test codomain(map_from_character_lattice_to_torusinvariant_weil_divisor_group(P2v2)) == torusinvariant_weil_divisor_group(P2v2)
        @test domain(map_from_cartier_divisor_group_to_picard_group(P2v2)) == cartier_divisor_group(P2v2)
        @test codomain(map_from_cartier_divisor_group_to_picard_group(P2v2)) == picard_group(P2v2)
        @test domain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(P2v2)) == cartier_divisor_group(P2v2)
        @test codomain(map_from_cartier_divisor_group_to_torusinvariant_divisor_group(P2v2)) == torusinvariant_weil_divisor_group(P2v2)
        @test coordinate_names(P2v2) == ["x1", "x2", "x3"]
    end

    @testset "Blowup" begin
        blowup_variety = blowup_on_ith_minimal_torus_orbit(P2, 1, "e")
        @test is_normal(blowup_variety) == true
        @test is_affine(blowup_variety) == false
        @test is_projective(blowup_variety) == true
        @test is_smooth(blowup_variety) == true
        @test is_complete(blowup_variety) == true
        @test has_torusfactor(blowup_variety) == false
        @test is_orbifold(blowup_variety) == true
        @test is_simplicial(blowup_variety) == true
        @test betti_number(blowup_variety, 0) == 1
        @test betti_number(blowup_variety, 1) == 0
        @test betti_number(blowup_variety, 2) == 2
        @test betti_number(blowup_variety, 3) == 0
        @test betti_number(blowup_variety, 4) == 1
        @test euler_characteristic(blowup_variety) == 4
        @test rank(picard_group(blowup_variety)) == 2
        @test is_projective_space(blowup_variety) == false
    end
end
