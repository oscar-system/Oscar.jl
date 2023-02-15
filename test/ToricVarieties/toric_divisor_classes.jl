using Oscar
using Test

@testset "Torus-invariant divisor classes (set_attributes = $set_attributes)" for set_attributes in [true, false]

    F5 = hirzebruch_surface(5; set_attributes)
    dP3 = del_pezzo_surface(3; set_attributes)
    P2 = projective_space(NormalToricVariety, 2; set_attributes)
    
    DC = ToricDivisorClass(F5, [fmpz(0), fmpz(0)])
    DC2 = ToricDivisorClass(F5, [1, 2])
    DC3 = ToricDivisorClass(dP3, [4, 3, 2, 1])
    DC4 = canonical_divisor_class(dP3)
    DC5 = anticanonical_divisor_class(dP3)
    DC6 = trivial_divisor_class(dP3)
    DC7 = ToricDivisorClass(P2, [1])
    DC8 = ToricDivisorClass(P2, [-1])
    
    @testset "Basic properties" begin
        @test is_trivial(toric_divisor(DC2)) == false
        @test is_effective(DC7) == true
        @test is_effective(DC8) == false
    end
    
    @testset "Basic attributes" begin
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
