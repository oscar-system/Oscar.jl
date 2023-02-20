using Oscar
using Test

@testset "Line bundles (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    dP1 = del_pezzo_surface(1; set_attributes)
    dP3 = del_pezzo_surface(3; set_attributes)
    
    l = ToricLineBundle(dP3, [1, 2, 3, 4])
    l2 = canonical_bundle(dP3)
    l3 = anticanonical_bundle(dP3)
    l4 = ToricLineBundle(dP3, trivial_divisor(dP3))
    l5 = ToricLineBundle(dP1, [fmpz(1), fmpz(2)])
    
    @testset "Should fail" begin
        @test_throws ArgumentError l * l5
    end
    
    @testset "Basic properties" begin
        @test is_trivial(l) == false
        @test is_basepoint_free(l) == false
        @test is_ample(l) == false
        @test is_very_ample(l) == false
    end
    
    @testset "Basic attributes" begin
        if set_attributes
            @test degree(l) == -6
            @test degree(l^(-1)) == 6
            @test degree(l*l) == -12
        else
            @test degree(l) == 10
            @test degree(l^(-1)) == -10
            @test degree(l*l) == 20
        end
        @test divisor_class(l).coeff == AbstractAlgebra.matrix(ZZ, [1 2 3 4])
        @test dim(toric_variety(l)) == 2
    end
    
    @testset "Arithmetic" begin
        @test (l == l5) == false
        @test (l == l2) == false
        @test (l2 * l3 == structure_sheaf(dP3)) == true
        @test (l * l4 * inv(l) == structure_sheaf(dP3)) == true
    end
end
