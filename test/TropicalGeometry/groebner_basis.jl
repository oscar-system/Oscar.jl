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

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector) - principal ideals" begin
        R,(x,y) = QQ["x","y"]
        I = ideal(R,[x^1+8*y^2])
        w = [-1,-1]
        nu = tropical_semiring_map(QQ)
        @test issetequal(groebner_basis(I,nu,w),gens(I))
    end

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector) - binomial ideals" begin
        R,(x,y) = QQ["x","y"]
        I = ideal(R,[x^1+8*y^2,x^2+y^4])
        w = [-2,-1]
        nu = tropical_semiring_map(QQ)
        @test issetequal(groebner_basis(I,nu,w),gens(I))
    end

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector) - linear ideals, trivial valuation" begin
        R,(x,y,z) = QQ["x","y","z"]
        nu = tropical_semiring_map(QQ)
        I = ideal(R,[x-z-1,y+z-1])
        wxyz = [1,2,3]
        @test issetequal(groebner_basis(I,nu,wxyz),[-y-z+1,x-y-2*z])

        wzyx = [3,2,1]
        groebner_basis(I,nu,wzyx)
        wzyx = [3,2,1]
        @test issetequal(groebner_basis(I,nu,wzyx),[-1//2*x-1//2*y+1,-1//2*x+1//2*y+z])
    end

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector) - trivial valuation, global ordering" begin

        R,(x,y) = QQ["x","y"]
        I = ideal(R,[x^1+8*y^2,x^2+y^4])
        w = [-1,-1]
        nu = tropical_semiring_map(QQ)
        groebner_basis(I,nu,w)
        @test issetequal(groebner_basis(I,nu,w),[8*y^2+x,x^2])
    end

end
