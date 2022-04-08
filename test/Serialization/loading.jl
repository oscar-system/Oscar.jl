@testset "loading" begin
    
    @testset "loading Vector{LinearProgram}" begin
        c = cube(3)
        LP0 = LinearProgram(c, [2,2,-3])
        LP1 = LinearProgram(c, [2,2,4])
        v = [LP0, LP1]
        loaded = load(joinpath(@__DIR__,"vlp.json"))
        @test length(v) == length(loaded)
        @test feasible_region(loaded[1]) == feasible_region(loaded[2])
        @test feasible_region(loaded[1]) == feasible_region(LP0)
        @test objective_function(loaded[1]) == objective_function(v[1])
        @test objective_function(loaded[2]) == objective_function(v[2])
        @test optimal_value(loaded[1]) == optimal_value(v[1])
        @test optimal_value(loaded[2]) == optimal_value(v[2])
    end
end
