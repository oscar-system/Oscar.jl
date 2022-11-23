@testset "Normal toric varieties" begin
    ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
    set_coordinate_names(ntv, ["x1", "x2", "y1", "y2"])
    ntv2 = NormalToricVariety(Oscar.cube(2))
    ntv3 = NormalToricVarietyFromGLSM(matrix(ZZ, [[1, 1, 1]]))
    ntv4 = NormalToricVarietiesFromStarTriangulations(convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]))
    ntv5 = NormalToricVariety(polarize(Polyhedron(Polymake.polytope.rand_sphere(5, 60; seed=42))))

    @testset "Basic properties" begin
        @test is_complete(ntv) == true
        @test rank(torusinvariant_cartier_divisor_group(ntv)) == 4
        @test rank(domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(ntv))) == 4
        @test is_complete(ntv2) == true
        @test length(ntv3) == 1
        @test is_projective_space(ntv3[1]) == true
        @test length(ntv4) == 2
        @test is_projective_space(ntv) == false
        @test is_projective_space(ntv2) == false
        @test is_projective_space(ntv3[1]) == true
        @test is_projective_space(ntv4[1]) == false
        @test is_projective_space(ntv4[2]) == false
        @test is_projective_space(ntv5) == false
    end

    @testset "Speed test for Stanley-Reisner ideal (at most a few seconds)" begin
        duration = @elapsed stanley_reisner_ideal(ntv5)
        @test duration < 10
        @test ngens(stanley_reisner_ideal(ntv5)) == 1648
    end
end
