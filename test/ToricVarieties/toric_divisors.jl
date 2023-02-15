using Oscar
using Test

@testset "Torus-invariant divisors (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    F5 = hirzebruch_surface(5; set_attributes)
    dP3 = del_pezzo_surface(3; set_attributes)
    P2 = projective_space(NormalToricVariety, 2; set_attributes)
    
    D=ToricDivisor(F5, [0, 0, 0, 0])
    D2 = DivisorOfCharacter(F5, [1, 2])
    D3 = ToricDivisor(dP3, [1, 0, 0, 0, 0, 0])
    D4 = canonical_divisor(dP3)
    D5 = anticanonical_divisor(dP3)
    D6 = trivial_divisor(dP3)
    D7 = ToricDivisor(P2, [1,1,0])
    D8 = ToricDivisor(P2, [1,-1,0])
    D9 = ToricDivisor(P2, [0,-1,0])
    
    @testset "Should fail" begin
        @test_throws ArgumentError ToricDivisor(F5, [0, 0, 0])
        @test_throws ArgumentError D+D3
        @test_throws ArgumentError D-D3
    end
    
    @testset "Basic properties" begin
        @test is_prime(D) == false
        @test is_prime(D2) == false
        @test is_prime(D3) == true
        @test is_cartier(D) == true
        @test is_cartier(D2) == true
        @test is_principal(D) == true
        @test is_principal(D2) == true
        @test is_trivial(D) == true
        @test is_trivial(D2) == false
        @test is_basepoint_free(D) == true
        @test is_basepoint_free(D2) == true
        @test is_ample(D) == false
        @test is_ample(D2) == false
        @test is_very_ample(D) == false
        @test is_very_ample(D2) == false
        @test is_nef(D) == true
        @test is_nef(D2) == true
        @test is_integral(D) == true
        @test is_integral(D2) == true
        @test is_q_cartier(D) == true
        @test is_q_cartier(D2) == true
        @test is_prime(D) == false
        @test is_prime(D2) == false
        @test is_effective(D7) ==  true
        @test is_effective(D8) ==  false
    end
    
    @testset "Basic attributes" begin
        @test coefficients(D) == [0, 0, 0, 0]
        @test coefficients(D2) == [1, 2, 9, -2]
        @test dim(toric_variety(D)) == 2
        @test dim(polyhedron(D)) == 0
        @test ambient_dim(polyhedron(D)) == 2
    end
    
    @testset "Arithmetic" begin
        @test (D == D2) == false
        @test (D4 + D5 == D6) == true
        @test is_principal(fmpz(2)*D+D2) == true
        @test is_principal(2*D-D2) == true
        @test coefficients(D2+D2) == coefficients(2*D2)
        @test coefficients(D2-D2) == [0, 0, 0, 0]
    end
end
