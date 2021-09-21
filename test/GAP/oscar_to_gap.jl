@testset "fmpz" begin
    # small (GAP) integer
    x = fmpz(17)
    val = 17
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # large GAP integer
    x = fmpz(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val
end

@testset "fmpq" begin
    # small (GAP) integer
    x = fmpz(17)
    val = 17
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # large GAP integer
    x = fmpz(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # non-integer rational, small numerator and denominator
    x = fmpq(2, 3)
    val = GAP.evalstr("2/3")
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # non-integer rational, large numerator and denominator
    x = fmpq(fmpz(2)^65, fmpz(3)^40)
    val = GAP.evalstr("2^65/3^40")
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val
end

@testset "fmpz_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.ZZ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # matrix containing small and large integers
    x = Nemo.ZZ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val
end

@testset "fmpq_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.QQ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # matrix containing small and large integers
    x = Nemo.QQ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # matrix containing non-integer rationals, small numerator and denominator
    x = Nemo.QQ[fmpq(1, 2) 2; 3 4]
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val

    # matrix containing non-integer rationals, large numerator and denominator
    x = Nemo.QQ[fmpq(fmpz(2)^65, fmpz(3)^40) 2; 3 4]
    val = GAP.evalstr( "[ [ 2^65/3^40, 2 ], [ 3, 4 ] ]" )
    @test GAP.julia_to_gap(x) == val
    @test GAP.GapObj(x) == val
end

@testset "GapGroup and GapGroupElem" begin
    # `GapGroup` to GAP group, Perm
    G = symmetric_group(5)
    val = GAP.evalstr("SymmetricGroup(5)")
    @test GAP.GapObj(G) == val

    # `GapGroupElem` to GAP group element, Perm
    g = perm(G, [2,3,1,5,4])
    val = GAP.evalstr("(1,2,3)(4,5)")
    @test GAP.GapObj(g) == val
end
