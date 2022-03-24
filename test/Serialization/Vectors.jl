@testset "Vectors" begin
    
    mktempdir() do path
        @testset "Vector{LinearProgram}" begin
            c = cube(3)
            LP0 = LinearProgram(c, [2,2,-3])
            LP1 = LinearProgram(c, [2,2,4])
            v = [LP0, LP1]
            save(v, joinpath(path, "vlp.json"))
            loaded = load(joinpath(path, "vlp.json"))
            @test length(v) == length(loaded)
            @test feasible_region(loaded[1]) == feasible_region(loaded[2])
            @test feasible_region(loaded[1]) == feasible_region(LP0)
            @test objective_function(loaded[1]) == objective_function(v[1])
            @test objective_function(loaded[2]) == objective_function(v[2])
            @test optimal_value(loaded[1]) == optimal_value(v[1])
            @test optimal_value(loaded[2]) == optimal_value(v[2])
        end

        @testset "Vector{gfp_fmpz_elem}" begin
            F = GF(fmpz(77777732222322222232222222223))
            one = F(1)
            minusone = F(-1)
            v = [one, minusone]
            save(v, joinpath(path, "vgfe.json"))
            loaded = load(joinpath(path, "vgfe.json"))
            @test v == loaded
        end
        
        @testset "Vector{gfp_elem}" begin
            F = GF(7)
            one = F(1)
            minusone = F(-1)
            v = [one, minusone]
            save(v, joinpath(path, "vge.json"))
            loaded = load(joinpath(path, "vge.json"))
            @test v == loaded
        end

        @testset "Vector{Any}" begin
            c = cube(3)
            LP0 = LinearProgram(c, [2,2,-3])
            v = [c, LP0]
            save(v, joinpath(path, "vany1.json"))
            loaded = load(joinpath(path, "vany1.json"))
            @test length(v) == length(loaded)
            @test loaded[1] isa Polyhedron
            @test loaded[2] isa LinearProgram
            @test loaded isa Vector{Any}
        end
        
        @testset "Vector{Union{Polyhedron, LinearProgram}}" begin
            c = cube(3)
            LP0 = LinearProgram(c, [2,2,-3])
            v = Vector{Union{Polyhedron, LinearProgram}}([c, LP0])
            save(v, joinpath(path, "vany2.json"))
            loaded = load(joinpath(path, "vany2.json"))
            @test length(v) == length(loaded)
            @test loaded[1] isa Polyhedron
            @test loaded[2] isa LinearProgram
            @test loaded isa Vector{Union{Polyhedron, LinearProgram}}
        end
    end
end
