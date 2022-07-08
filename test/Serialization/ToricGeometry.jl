@testset "ToricGeometry" begin
    
    mktempdir() do path
        @testset "NormalToricVariety" begin
            pp = projective_space(NormalToricVariety, 2)
            test_save_load_roundtrip(path, pp) do loaded
              @test rays(pp) == rays(loaded)
              @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
            end
        end

        @testset "ToricDivisor" begin
            pp = projective_space(NormalToricVariety, 2)
            td0 = ToricDivisor(pp, [1,1,2])
            td1 = ToricDivisor(pp, [1,1,3])
            vtd = [td0, td1]
            test_save_load_roundtrip(path, vtd) do loaded
              @test coefficients(td0) == coefficients(loaded[1])
              @test coefficients(td1) == coefficients(loaded[2])
              @test toric_variety(loaded[1]) == toric_variety(loaded[2])
            end
        end
    end

end
