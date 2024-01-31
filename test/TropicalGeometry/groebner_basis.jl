@testset "src/TropicalGeometry/groebner_basis.jl" begin

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector)" begin
        R,(x,y) = QQ["x","y"]
        I = ideal(R,[x^2+8*y^2,x+y])
        w = [1,-1]
        nuMin = tropical_semiring_map(QQ)
        nuMax = tropical_semiring_map(QQ,max)
        @test issetequal(groebner_basis(I,nuMin,w),[x^2,x+y]) # trivial, min
        @test issetequal(groebner_basis(I,nuMax,w),[y^2,x+y]) # trivial, max

        nuMin = tropical_semiring_map(QQ,2)
        nuMax = tropical_semiring_map(QQ,2,max)
        @test issetequal(groebner_basis(I,nuMin,w),[x^2-8*x*y,x+y]) # padic, min
        @test issetequal(groebner_basis(I,nuMax,w),[y^2,x+y]) # padic, max

        K,t = rational_function_field(GF(2),"t")
        R,(x,y) = K["x","y"]
        I = ideal(R,[x^2+t^3*y^2,x+y])
        nuMin = tropical_semiring_map(K,t)
        nuMax = tropical_semiring_map(K,t,max)
        @test issetequal(groebner_basis(I,nuMin,w),[x^2-t^3*x*y,x+y]) # tadic, min
        @test issetequal(groebner_basis(I,nuMax,w),[y^2,x+y]) # tadic, max
    end

end
