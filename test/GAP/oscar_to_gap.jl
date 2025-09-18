@testset "single cyclotomics" begin
    # from cyclotomic fields
    F, z = cyclotomic_field(5)
    e5 = GAP.Globals.E(5)
    @test GAP.Obj(z^2+z+1) == e5^2 + e5 + 1

    # from `QQAbFieldElem`
    F, z = abelian_closure(QQ)
    @test GAP.Obj(z(5) + z(5)^4) == e5 + e5^4

    # not supported conversions
    F, z = quadratic_field(5)
    @test_throws ArgumentError GAP.Obj(z)
end

@testset "matrices over a cyclotomic field" begin
    F, z = cyclotomic_field(5)
    @test GAP.Obj(matrix(F, 2, 2, [z^0, z, z^2, z^3])) == GAP.evalstr("[ [ 1, E(5) ], [ E(5)^2, E(5)^3 ] ]")

    F, z = quadratic_field(5)
    @test_throws ArgumentError GAP.Obj(matrix(F, 1, 1, [z]))
end

@testset "GapGroup and GapGroupElem" begin
    # `GapGroup` to GAP group, Perm
    G = symmetric_group(5)
    val = GAP.Globals.SymmetricGroup(5)
    @test GAP.Obj(G) == val

    # `GapGroupElem` to GAP group element, Perm
    g = perm(G, [2,3,1,5,4])
    val = GAP.evalstr("(1,2,3)(4,5)")
    @test GAP.Obj(g) == val
end
