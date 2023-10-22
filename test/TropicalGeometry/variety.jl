@testset "src/TropicalGeometry/{variety,hypersurface,curve,linear_space}.jl" begin

    # constructing various tropicalizations two ways
    # and comparing their results
    @testset "hypersurface and linear spaces" begin
        R,(x,y,z) = QQ["x","y","z"]
        f = x+2*y+4*z
        nu = tropical_semiring_map(QQ,2)
        TropH = tropical_hypersurface(f,nu)
        TropL = tropical_linear_space(ideal(R,f),nu)
        @test issetequal(maximal_polyhedra(TropH),maximal_polyhedra(TropL))
        nu = tropical_semiring_map(QQ,2,max)
        TropH = tropical_hypersurface(f,nu)
        TropL = tropical_linear_space(ideal(R,f),nu)
        @test issetequal(maximal_polyhedra(TropH),maximal_polyhedra(TropL))
    end

    @testset "binomial ideals and linear spaces" begin
        R,(x,y,z,w) = QQ["x","y","z","w"]
        I = ideal(R,[x+2*y,z+4*w])
        nu = tropical_semiring_map(QQ,2)
        TropV = first(tropical_variety(I,nu))
        TropL = tropical_linear_space(I,nu)
        @test issetequal(maximal_polyhedra(TropV),maximal_polyhedra(TropL))
        nu = tropical_semiring_map(QQ,2,max)
        TropV = first(tropical_variety(I,nu))
        TropL = tropical_linear_space(I,nu)
        @test issetequal(maximal_polyhedra(TropV),maximal_polyhedra(TropL))
    end

    @testset "hypersurface and binomial ideals" begin
        R,(x,y,z) = QQ["x","y","z"]
        f = x*y*z+2
        nu = tropical_semiring_map(QQ,2)
        TropH = tropical_hypersurface(f,nu)
        TropV = first(tropical_variety(ideal(R,f),nu))
        @test issetequal(maximal_polyhedra(TropH),maximal_polyhedra(TropV))
        nu = tropical_semiring_map(QQ,2,max)
        TropV = first(tropical_variety(ideal(R,f),nu))
        TropH = tropical_hypersurface(f,nu)
        @test issetequal(maximal_polyhedra(TropH),maximal_polyhedra(TropV))
    end

end
