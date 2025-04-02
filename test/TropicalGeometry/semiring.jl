@testset "src/TropicalGeometry/semiring.jl" begin

    @testset "constructors and convention" begin
        @test convention(tropical_semiring())==min
        @test convention(tropical_semiring(min))==min
        @test convention(tropical_semiring(max))==max
    end

    for minOrMax in (min,max)
        T = tropical_semiring(minOrMax)
        @testset "conversions" begin
            for K in (QQ,ZZ,Int)
                @test T(K(0)) == one(T)
                @test K(one(T)) == 0
                @test K(T(1)) == 1
                @test K(T(1),preserve_ordering=true) == (minOrMax==min ? 1 : -1)
            end
        end

        @testset "arithmetics" begin
            @test zero(T)+one(T) == one(T)
            @test zero(T)*one(T) == zero(T)
            @test T(0)+T(1) == (minOrMax==min ? T(0) : T(1))
        end

        @testset "polynomial exponentiation" begin
            Tx, x = polynomial_ring(T, 1)
            f = x[1]^4 + 3*x[1]
            @test f^3 == f*f*f
        end
    end

end
