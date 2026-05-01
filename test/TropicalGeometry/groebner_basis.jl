@testset "src/TropicalGeometry/groebner_basis.jl" begin

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::AbstractVector)" begin
        R,(x,y) = QQ[:x, :y]
        I = ideal(R,[x^2+x*y+8*y^2,x+y])
        w = [1,-1]
        nuMin = tropical_semiring_map(QQ)
        nuMax = tropical_semiring_map(QQ,max)
        @test issetequal(groebner_basis(I,nuMin,w),[x^2,x+y]) # trivial, min
        @test issetequal(groebner_basis(I,nuMax,w),[y^2,x+y]) # trivial, max

        @test issetequal(groebner_basis(I,nuMin,point_vector(w)),[x^2,x+y]) # testing PointVector
        @test issetequal(groebner_basis(I,nuMin,ray_vector(w)),[x^2,x+y]) # testing RayVector


        nuMin = tropical_semiring_map(QQ,2)
        nuMax = tropical_semiring_map(QQ,2,max)
        @test issetequal(groebner_basis(I,nuMin,w),[x^2,x+y]) # padic, min
        @test issetequal(groebner_basis(I,nuMax,w),[y^2,x+y]) # padic, max

        K,t = rational_function_field(GF(2),"t")
        R,(x,y) = K[:x, :y]
        I = ideal(R,[x^2+x*y+t^3*y^2,x+y])
        nuMin = tropical_semiring_map(K,t)
        nuMax = tropical_semiring_map(K,t,max)
        @test issetequal(groebner_basis(I,nuMin,w),[x^2,x+y]) # tadic, min
        @test issetequal(groebner_basis(I,nuMax,w),[y^2,x+y]) # tadic, max
    end

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::AbstractVector) - principal ideals" begin
        R,(x,y) = QQ[:x, :y]
        I = ideal(R,[x^1+8*y^2])
        w = [-1,-1]
        nu = tropical_semiring_map(QQ)
        @test issetequal(groebner_basis(I,nu,w),gens(I))
    end

    @testset "groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::AbstractVector) - weight adjustments for homogeneous ideals" begin
        R, (x1,x2,x3) = polynomial_ring(QQ,3)
        nu = tropical_semiring_map(QQ,max)
        nu2 = tropical_semiring_map(QQ,2,max)
        I = ideal(2*x1+x2+x3, 2*x1+3*x2+4*x3)
        w = [1,1,typemax(Int16)+1]
        @test begin
            # entries of w are too large for Singular but can be made to fit Singular for homogeneous ideals by translating w.r.t. all ones vector
            # this checks whether this is done by testing whether the groebner_basis calls below run without raising an error
            groebner_basis(I,nu,w)
            groebner_basis(I,nu2,w)
        end
    end

end
