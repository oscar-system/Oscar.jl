@testset "src/TropicalGeometry/semiring.jl" begin

    @testset "constructors and convention" begin
        @test convention(tropical_semiring())==Min
        @test convention(tropical_semiring(Min))==Min
        @test convention(tropical_semiring(Max))==Max
    end

    for MinOrMax in (Min,Max)
        T = tropical_semiring(MinOrMax)
        @testset "conversions" begin
            for K in (QQ,ZZ,Int)
                @test T(K(0)) == one(T)
                @test K(one(T)) == 0
                @test K(T(1)) == 1
                @test K(T(1),preserve_ordering=true) == (MinOrMax==Min ? 1 : -1)
            end
        end

        @testset "arithmetics" begin
            @test zero(T)+one(T) == one(T)
            @test zero(T)*one(T) == zero(T)
            @test T(0)+T(1) == (MinOrMax==Min ? T(0) : T(1))
        end
    end

end
