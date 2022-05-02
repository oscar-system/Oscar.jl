@testset "GAP matrix groups" begin

    m = matrix(ZZ, [0 1 ; -1 0])
    gapM = wrap_for_gap(m)
    G = GAP.Globals.Group(gapM)
    GAP.Globals.SetNiceMorphismForJuliaMatrixRepGroup(G)
    @test GAP.Globals.Size(G) == 4
    
end
