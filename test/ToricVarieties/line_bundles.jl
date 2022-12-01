@testset "Toric line bundles" begin
    dP3 = NormalToricVariety([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]])
    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    D2 = DivisorOfCharacter(F5, [1, 2])
    dP1 = NormalToricVariety([[1, 0], [0, 1], [-1, 0], [-1, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])

    l = ToricLineBundle(dP3, [1, 2, 3, 4])
    l2 = ToricLineBundle(D2)
    l3 = ToricLineBundle(dP1, [fmpz(1), fmpz(2)])
    l4 = canonical_bundle(dP3)
    l5 = anticanonical_bundle(dP3)
    l6 = ToricLineBundle(dP3, trivial_divisor(dP3))
    l7 = structure_sheaf(dP3)

    @testset "Should fail" begin
        @test_throws ArgumentError l * l3
    end

    @testset "Basic properties" begin
        @test is_trivial(l) == false
        @test is_basepoint_free(l) == false
        @test is_ample(l) == false
        @test is_very_ample(l) == false
        @test degree(l) == 10
        @test divisor_class(l).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
        @test dim(toric_variety(l)) == 2
    end

    @testset "Arithmetic" begin
        @test (l == l3) == false
        @test (l4 * l5 == l7) == true
        @test (l * l6 * inv(l) == l7) == true
        @test degree(l^(-1)) == -10
        @test degree(l*l) == 20
    end
end
