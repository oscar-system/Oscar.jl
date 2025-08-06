@testset "src/TropicalGeometry/intersection.jl" begin

    # stably intersecting two tropical hypersurfaces
    @testset "hypersurface and linear spaces" begin
        R,(x,y,z) = QQ[:x, :y, :z]
        f = x+2*y+4*z
        TropH = tropical_hypersurface(f)
        TropL = tropical_linear_space(ideal(R,f))
        TropHH = stable_intersection(TropH,TropH)
        TropLL = stable_intersection(TropL,TropL)
        @test issetequal(maximal_polyhedra(TropHH),maximal_polyhedra(TropLL))

        TropHHH = stable_intersection(TropHH,TropH)
        @test dim(TropHHH)<0

        TropLLL = stable_intersection(TropLL,TropL)
        @test dim(TropLLL)<0
    end

end
