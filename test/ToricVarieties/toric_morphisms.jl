using Oscar
using Test

@testset "Toric morphisms (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    source = projective_space(NormalToricVariety, 1; set_attributes)
    target = hirzebruch_surface(2; set_attributes)
    tm1 = ToricMorphism(source, matrix(ZZ, [[0, 1]]), target)
    tm2 = ToricIdentityMorphism(target)
    tm3 = ToricIdentityMorphism(hirzebruch_surface(4; set_attributes))
    
    @testset "Argument errors for toric morphisms" begin
        @test_throws ArgumentError tm1 * tm1
    end
    
    @testset "Arithmetics of toric morphisms" begin
        @test (tm1 == tm2) == false
        @test (tm2 + tm2 == 2*tm2) == true
        @test (tm2 + tm2 == fmpz(2)*tm2) == true
        @test (tm2 * tm2 == tm2) == true
    end
    
    @testset "Basic attributes of toric morphisms" begin
        @test (image(tm1) == codomain(tm1)) == false
        @test codomain(tm1) == domain(tm2)
        @test matrix(grid_morphism(tm2)) == matrix(ZZ, [[1, 0], [0, 1]])
        @test morphism_on_torusinvariant_weil_divisor_group(tm2).map == identity_matrix(ZZ, 4)
        @test morphism_on_torusinvariant_cartier_divisor_group(tm2).map == identity_matrix(ZZ, 4)
        @test grid_morphism(morphism_from_cox_variety(source)).map == matrix(ZZ, [[1], [-1]])
        @test is_affine(cox_variety(source)) == false
        @test matrix(morphism_on_class_group(tm3)) == matrix(ZZ, [[1, 0], [0, 1]])
        @test matrix(morphism_on_picard_group(tm3)) == matrix(ZZ, [[1, 0], [0, 1]])
    end
end
