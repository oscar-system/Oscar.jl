@testset "del Pezzo surfaces" begin

    dP0 = NormalToricVariety(normal_fan(Oscar.simplex(2)))
    dP1 = NormalToricVariety([[1, 0], [0, 1], [-1, 0], [-1, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    dP2 = NormalToricVariety([[1, 0], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]])
    dP3 = NormalToricVariety([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]])
    set_coordinate_names(dP3, ["x1", "e1", "x2", "e3", "x3", "e2"])
    dP3v2 = del_pezzo_surface(3)

    @testset "Should fail" begin
        @test_throws ArgumentError del_pezzo_surface(-1)
        @test_throws ArgumentError del_pezzo_surface(4)
    end

    @testset "Basic properties" begin
        @test length(torusinvariant_prime_divisors(dP0)) == 3
        @test length(torusinvariant_prime_divisors(dP1)) == 4
        @test length(torusinvariant_prime_divisors(dP2)) == 5
        @test rank(torusinvariant_cartier_divisor_group(dP3)) == 6
        @test length(torusinvariant_prime_divisors(dP3)) == 6
        @test rank(picard_group(dP3)) == 4
        @test picard_group(dP3) == codomain(map_from_torusinvariant_cartier_divisor_group_to_picard_group(dP3))
        @test is_projective_space(dP0) == true
        @test is_projective_space(dP1) == false
        @test is_projective_space(dP2) == false
        @test is_projective_space(dP3) == false
    end

    @testset "Constructors" begin
      phi = map_from_torusinvariant_weil_divisor_group_to_class_group(dP3v2)
      @test vcat([phi(x).coeff for x in gens(torusinvariant_weil_divisor_group(dP3v2))]) == matrix(ZZ, [[1, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
      @test coordinate_names(dP3v2) == ["x1", "x2", "x3", "e1", "e2", "e3"]
    end
end
