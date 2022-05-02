@testset "GAP matrix groups" begin

    m = matrix(ZZ, [0 1 ; -1 0])
    gapM = Oscar.MatrixGroups._wrap_for_gap(m)
    G = GAP.Globals.Group(gapM)
    GAP.Globals.SetNiceMorphismForJuliaMatrixRepGroup(G)
    @test GAP.Globals.Size(G) == 4
    
    m2 = matrix(QQ, [ -1 0; 0 1]);
    G = Oscar.MatrixGroups.MatrixGroup([m,m2])
    @test GAP.Globals.Size(G) == 8
end
