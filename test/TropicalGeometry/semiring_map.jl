@testset "src/TropicalGeometry/semiringMaps.jl" begin

    ###
    # constructing all possible semiring maps
    # and perforMing basic sanity checks
    ###
    for MinOrMax in (Min,Max)
        T = tropical_semiring(MinOrMax)

        @testset "trivial valuation" begin
            nu = tropical_semiring_map(GF(2),MinOrMax)
            @test nu(0)==zero(T)
            @test nu(1)==one(T)
        end

        @testset "p-adic valuation" begin
            nu = tropical_semiring_map(QQ,2,MinOrMax)
            @test nu(0)==zero(T)
            @test nu(1+2)==one(T)
            @test nu(2+4)==(MinOrMax==Min ? T(1) : T(-1))
        end

        @testset "t-adic valuation" begin
            K,t = rational_function_field(GF(2),"t")
            nu = tropical_semiring_map(K,t,MinOrMax)
            @test nu(0)==zero(T)
            @test nu(1+t)==one(T)
            @test nu(t+t^2)==(MinOrMax==Min ? T(1) : T(-1))
        end
    end

end
