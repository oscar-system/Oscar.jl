using Oscar
using Test
include("../src/WordCalculations.jl")

@testset "Test CalculateBetas for A2" begin
    type = "A"
    rank = 2
    word = [1, 2, 1]
    lie_algebra = BasisLieHighestWeight.LieAlgebraStructure(type, rank)
    
    betas = compute_betas(lie_algebra, word)
    
    # Expected beta values
    expected_betas = [
        [2, -1],
        [1, 1],
        [-1, 2]
    ]
    @test betas == expected_betas
end
