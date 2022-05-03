@testset "GAP matrix groups" begin

    m = matrix(ZZ, [0 1 ; -1 0])
    gapM = Oscar.MatrixGroups._wrap_for_gap(m)
    G = GAP.Globals.Group(gapM)
    GAP.Globals.SetNiceMorphismForJuliaMatrixRepGroup(G)
    @test GAP.Globals.Size(G) == 4
    
    m2 = matrix(QQ, [ -1 0; 0 1])
    G = Oscar.MatrixGroups.MatrixGroup([m,m2])
    @test_throws ERROR: MethodError: no method matching MatrixGroup(::Vector{MatElem})
    
    m1 = matrix(QQ, [0 1 ; -1 0])
    G = Oscar.MatrixGroups.MatrixGroup([m1,m2])
    @test GAP.Globals.Size(G) == 8
    
    K, a = quadratic_field(-1)
    m1 = matrix(K, [ a 0 ; 0 a ])
    m2 = matrix(K, [ 0 1 ; -1 0 ])
    G = Oscar.MatrixGroups.MatrixGroup([m1,m2])
    @test GAP.Globals.Size(G) == 8
end
