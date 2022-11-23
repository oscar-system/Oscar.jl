@testset "Torus invariant divisor classes" begin

    F5 = NormalToricVariety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    dP3 = NormalToricVariety([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]])
    DC = ToricDivisorClass(F5, [fmpz(0), fmpz(0)])
    DC2 = ToricDivisorClass(F5, [1, 2])
    DC3 = ToricDivisorClass(dP3, [4, 3, 2, 1])
    DC4 = canonical_divisor_class(dP3)
    DC5 = anticanonical_divisor_class(dP3)
    DC6 = trivial_divisor_class(dP3)

    @testset "Basic properties" begin
        @test is_trivial(toric_divisor(DC2)) == false
        @test rank(parent(divisor_class(DC2))) == 2
        @test dim(toric_variety(DC2)) == 2
    end

    @testset "Arithmetic" begin
        @test is_trivial(fmpz(2)*DC+DC2) == false
        @test is_trivial(2*DC-DC2) == false
        @test (DC == DC2) == false
        @test (DC4 - DC5 == DC6) == false
        @test (DC == DC3) == false
    end
end
