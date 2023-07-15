using Oscar
using Test

include("../src/VectorSpaceBases.jl")

@testset "Test VectorSpaceBases" begin
    @testset "1" begin
        @test isequal(1, 1)
    end
end
