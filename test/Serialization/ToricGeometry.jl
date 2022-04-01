@testset "ToricGeometry" begin
    
    mktempdir() do path
        @testset "NormalToricVariety" begin
            pp = projective_space(NormalToricVariety, 2)
            filename = joinpath(path, "pp.ntv")
            save(pp, filename)
            loaded = load(filename)
            @test rays(pp) == rays(loaded)
            @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
        end
    end
end
