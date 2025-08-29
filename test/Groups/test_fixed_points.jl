using Test
using Oscar

@testset "fixed_points" begin
    # Full symmetric group -- no points fixed
    G = symmetric_group(5)
    @test fixed_points(G) == Int[]
end








