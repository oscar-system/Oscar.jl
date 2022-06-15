@testset "fmpz" begin
    # small (GAP) integer
    x = fmpz(17)
    val = 17
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # large GAP integer
    x = fmpz(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val
end

@testset "fmpq" begin
    # small (GAP) integer
    x = fmpz(17)
    val = 17
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # large GAP integer
    x = fmpz(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # non-integer rational, small numerator and denominator
    x = fmpq(2, 3)
    val = GAP.evalstr("2/3")
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # non-integer rational, large numerator and denominator
    x = fmpq(fmpz(2)^65, fmpz(3)^40)
    val = GAP.evalstr("2^65/3^40")
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val
end

@testset "fmpz_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.ZZ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # matrix containing small and large integers
    x = Nemo.ZZ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val
end

@testset "fmpq_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.QQ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # matrix containing small and large integers
    x = Nemo.QQ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # matrix containing non-integer rationals, small numerator and denominator
    x = Nemo.QQ[fmpq(1, 2) 2; 3 4]
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val

    # matrix containing non-integer rationals, large numerator and denominator
    x = Nemo.QQ[fmpq(fmpz(2)^65, fmpz(3)^40) 2; 3 4]
    val = GAP.evalstr( "[ [ 2^65/3^40, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.Obj(x) == val
end

@testset "single cyclotomics" begin
    # from cyclotomic fields
    F, z = CyclotomicField(5)
    e5 = GAP.Globals.E(5)
    @test GAP.Obj(z^2+z+1) == e5^2 + e5 + 1

    # from `QQAbElem`
    F, z = abelian_closure(QQ)
    @test GAP.Obj(z(5) + z(5)^4) == e5 + e5^4

    # not supported conversions
    F, z = quadratic_field(5)
    @test_throws ArgumentError GAP.Obj(z)
end

@testset "matrices over a cyclotomic field" begin
    F, z = CyclotomicField(5)
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

@testset "Set($coll)" for coll in [
                            [7, 1, 5, 3, 10],
                            fmpz[7, 1, 5, 3, 10],
                            [:c,:b,:a,:b],
                            [ (1,:a), (1,:b), (2,:a), (2,:b) ],
                        ]
    # create a set
    s = Set(coll)
    # create sort duplicate free list (this is how GAP represents sets)
    l = sort(unique(coll))

    x = GAP.GapObj(s)
    @test x == GAP.GapObj(l)

    x = GAP.GapObj(s; recursive=true)
    @test x == GAP.GapObj(l; recursive=true)
end
