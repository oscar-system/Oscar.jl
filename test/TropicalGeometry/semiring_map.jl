@testset "src/TropicalGeometry/semiringMaps.jl" begin

    ###
    # constructing all possible semiring maps
    # and performing basic sanity checks
    ###
    for minOrMax in (min,max)
        T = tropical_semiring(minOrMax)

        @testset "trivial valuation" begin
            nu = tropical_semiring_map(GF(2),minOrMax)
            @test nu(0)==zero(T)
            @test nu(1)==one(T)
        end

        @testset "p-adic valuation" begin
            nu = tropical_semiring_map(QQ,2,minOrMax)
            @test nu(0)==zero(T)
            @test nu(1+2)==one(T)
            @test nu(2+4)==(minOrMax==min ? T(1) : T(-1))
        end

        @testset "t-adic valuation" begin
            K,t = rational_function_field(GF(2),"t")
            nu = tropical_semiring_map(K,t,minOrMax)
            @test nu(0)==zero(T)
            @test nu(1+t)==one(T)
            @test nu(t+t^2)==(minOrMax==min ? T(1) : T(-1))
        end
    end

end
