@testset "ToricGeometry" begin
    
    mktempdir() do path
        @testset "NormalToricVariety" begin
            pp = projective_space(NormalToricVariety, 2)
            filename = joinpath(path, "pp.ntv")
            save(filename, pp)
            loaded = load(filename)
            @test rays(pp) == rays(loaded)
            @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
        end

        @testset "ToricDivisor" begin
            pp = projective_space(NormalToricVariety, 2)
            filename = joinpath(path, "td.json")
            td0 = ToricDivisor(pp, [1,1,2])
            td1 = ToricDivisor(pp, [1,1,3])
            vtd = [td0, td1]
            save(filename, vtd)
            loaded = load(filename)
            @test coefficients(td0) == coefficients(loaded[1])
            @test coefficients(td1) == coefficients(loaded[2])
            @test toric_variety(loaded[1]) == toric_variety(loaded[2])
        end
    end

end
