@testset "GAP matrix groups" begin
    m = matrix(ZZ, [0 1 ; -1 0])
    gapM = Oscar.MatrixGroups._wrap_for_gap(m)
    G = GAP.Globals.Group(gapM)
    GAP.Globals.SetNiceMorphismForJuliaMatrixRepGroup(G)
    @test GAP.Globals.Size(G) == 4

    m2 = matrix(QQ, [ -1 0; 0 1])
    @test_throws MethodError Oscar.MatrixGroups.matrix_group([m, m2])

    m1 = matrix(QQ, [0 1 ; -1 0])
    G = Oscar.MatrixGroups.matrix_group([m1, m2])
    @test GAP.Globals.Size(G) == 8

    K, a = quadratic_field(-1)
    m1 = matrix(K, [a 0 ; 0 a])
    m2 = matrix(K, [0 1 ; -1 0])
    G = Oscar.MatrixGroups.matrix_group([m1, m2])
    @test GAP.Globals.Size(G) == 8
    Ggens = GAP.Globals.GeneratorsOfGroup(G)
    c = GAP.Globals.Centralizer(G, GAP.Globals.Product(Ggens))
    @test GAP.Globals.Size(c) == 8

    # A5 as a matrix group in characteristic zero
    m1 = matrix(QQ, [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; -1 -1 -1 -1])
    m2 = matrix(QQ, [0 1 0 0 ; 0 0 0 1 ; 0 0 1 0 ; 1 0 0 0])
    G = Oscar.MatrixGroups.matrix_group([m1, m2])
    @test GAP.Globals.Size(G) == 60
    Ggens = GAP.Globals.GeneratorsOfGroup(G)
    x = GAP.Globals.Product(Ggens)
    c = GAP.Globals.Centralizer(G, x)
    @test GAP.Globals.Size(c) == 5
    iso = GAP.Globals.NiceMonomorphism(G)
    y = GAP.Globals.ImagesRepresentative(iso, x)
    @test x == GAP.Globals.PreImagesRepresentative(iso, y)

    # infinite group
    @test_throws ErrorException Oscar.MatrixGroups.matrix_group(
                                  [matrix(ZZ, [1 1 ; 0 1])])
end
