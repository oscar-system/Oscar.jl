using Oscar
using Test

include("../src/WordCalculations.jl")

@testset "Test CalculateBetas for A2" begin
    type = "A"
    rank = 2
    word = [1, 2, 1]
    
    betas = compute_betas(type, rank, word)
    betas = [Vector{Int}(GAP.Globals.List(b)) for b in betas]
    
    # Expected beta values
    expected_betas = [
        [2, -1],
        [1, 1],
        [-1, 2]
    ]
    @test betas == expected_betas
end
