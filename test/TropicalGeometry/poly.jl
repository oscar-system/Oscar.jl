@testset "src/TropicalGeometry/poly.jl" begin
    @testset "tropical_polynomial(::MPolyRingElem,::TropicalSemiringMap)" begin
        K,t = rational_function_field(GF(2),"t")
        nu = tropical_semiring_map(K,t)
        R,(x,y) = K[:x, :y]
        f = x+t*y+t^2
        tropf = tropical_polynomial(f,nu)
        @test issetequal(coefficients(tropf),tropical_semiring(nu).([0,1,2]))
        @test issetequal(exponents(tropf),[[1,0],[0,1],[0,0]])
    end
end
