@testset "src/TropicalGeometry/initial.jl" begin

    @testset "initial(::MPolyRingElem,::TropicalSemiringMap,::Vector)" begin
        R,(x,y) = QQ["x","y"]
        f = x^2+y^2+2*x+2*y
        w = [-1,-1]
        nuMin = tropical_semiring_map(QQ)
        nuMax = tropical_semiring_map(QQ,max)
        @test initial(f,nuMin,w) == x^2+y^2 # trivial, min
        @test initial(f,nuMax,w) == 2*x+2*y # trivial, max

        nuMin = tropical_semiring_map(QQ,2)
        nuMax = tropical_semiring_map(QQ,2,max)
        S,(x,y) = GF(2)["x","y"]
        @test initial(f,nuMin,w) == x^2+y^2 # padic, min
        @test initial(f,nuMax,w) == x^2+y^2+x+y # tadic, max

        K,t = rational_function_field(GF(2),"t")
        R,(x,y) = K["x","y"]
        f = x^2+y^2+t*x+t*y
        w = [-1,-1]
        nuMin = tropical_semiring_map(K,t)
        nuMax = tropical_semiring_map(K,t,max)
        S,(x,y) = GF(2)["x","y"]
        @test initial(f,nuMin,w) == x^2+y^2 # tadic, min
        @test initial(f,nuMax,w) == x^2+y^2+x+y # tadic, max
    end

    # no tests for initial(::MPolyIdeal,::TropicalSemiringMap,::Vector),
    # function is a simple composition of
    # - initial(::MPolyRingElem,::TropicalSemiringMap,::Vector) and
    # - groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector)

end
