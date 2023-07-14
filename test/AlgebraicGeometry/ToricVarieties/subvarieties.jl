using Oscar
using Test

@testset "Closed subvarieties (set_attributes = $set_attributes)" for set_attributes in [true, false]
    
    antv = normal_toric_variety([[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]], [[1, 2, 3, 4]])
    
    ntv = normal_toric_variety(cube(2); set_attributes)
    (x1, x2, y1, y2) = gens(cox_ring(ntv));
    sv1 = closed_subvariety_of_toric_variety(ntv, [x1])
    sv2 = closed_subvariety_of_toric_variety(ntv, [x1^2+x1*x2+x2^2, y2])
    
    P3 = projective_space(NormalToricVariety, 3; set_attributes)
    sv3 = closed_subvariety_of_toric_variety(P3, [gens(cox_ring(P3))[1]^2])
    
    @testset "Should fail" begin
        @test_throws ArgumentError closed_subvariety_of_toric_variety(ntv, [x1 - y1])
        @test_throws ArgumentError closed_subvariety_of_toric_variety(antv, [gens(cox_ring(antv))[1]])
    end
    
    @testset "Basic properties" begin
        @test is_empty(sv1) == false
        @test is_empty(sv2) == false
        @test is_empty(sv3) == false
        @test radical(sv1) == defining_ideal(sv1)
        @test dim(toric_variety(sv1)) == 2
    end
end
