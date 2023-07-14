using Oscar
using Test

@testset "Direct products" begin
    
    F5 = normal_toric_variety([[1, 0], [0, 1], [-1, 5], [0, -1]], [[1, 2], [2, 3], [3, 4], [4, 1]])
    P2 = normal_toric_variety(normal_fan(Oscar.simplex(2)))
    variety = F5 * P2
    
    @testset "Properties of a direct product of two toric varieties" begin
        @test is_smooth(variety) == true
        @test is_complete(variety) == true
        @test has_torusfactor(variety) == false
    end
        
    @testset "Betti numbers of a direct product of two toric varieties" begin
        @test betti_number(variety, 0) == 1
        @test betti_number(variety, 1) == 0
        @test betti_number(variety, 2) == 3
        @test betti_number(variety, 3) == 0
        @test betti_number(variety, 4) == 4
        @test betti_number(variety, 5) == 0
        @test betti_number(variety, 6) == 3
        @test betti_number(variety, 7) == 0
        @test betti_number(variety, 8) == 1
    end
end
