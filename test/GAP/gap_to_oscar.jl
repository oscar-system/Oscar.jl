@testset "fmpz" begin
    # small (GAP) integer
    x = fmpz(17)
    val = 17
    @test GAP.gap_to_julia(fmpz, val) == x
    @test convert(fmpz, val) == x
    @test fmpz(val) == x
    @test ZZ(val) == x

    # large GAP integer
    x = fmpz(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.gap_to_julia(fmpz, val) == x
    @test convert(fmpz, val) == x
    @test fmpz(val) == x
    @test ZZ(val) == x

    # non-integer
    val = GAP.evalstr("1/2")
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpz, val)
end

@testset "fmpq" begin
    # small (GAP) integer
    x = fmpq(17)
    val = 17
    @test GAP.gap_to_julia(fmpq, val) == x
    @test convert(fmpq, val) == x
    @test fmpq(val) == x
    @test QQ(val) == x

    # large GAP integer
    x = fmpq(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.gap_to_julia(fmpq, val) == x
    @test convert(fmpq, val) == x
    @test fmpq(val) == x
    @test QQ(val) == x

    # non-integer rational, small numerator and denominator
    x = fmpq(2, 3)
    val = GAP.evalstr("2/3")
    @test GAP.gap_to_julia(fmpq, val) == x
    @test convert(fmpq, val) == x
    @test fmpq(val) == x
    @test QQ(val) == x

    # non-integer rational, large numerator and denominator
    x = fmpq(fmpz(2)^65, fmpz(3)^40)
    val = GAP.evalstr("2^65/3^40")
    @test GAP.gap_to_julia(fmpq, val) == x
    @test convert(fmpq, val) == x
    @test fmpq(val) == x
    @test QQ(val) == x

    # non-rational
    val = GAP.evalstr("()")
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpq, val)
end

@testset "fmpz_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.ZZ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpz_mat, val) == x
    @test convert(fmpz_mat, val) == x
    @test fmpz_mat(val) == x

    # matrix containing small and large integers
    x = Nemo.ZZ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpz_mat, val) == x
    @test convert(fmpz_mat, val) == x
    @test fmpz_mat(val) == x

    # matrix containing non-integers
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpz_mat, val)
end

@testset "fmpq_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.QQ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test convert(fmpq_mat, val) == x
    @test fmpq_mat(val) == x

    # matrix containing small and large integers
    x = Nemo.QQ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test convert(fmpq_mat, val) == x
    @test fmpq_mat(val) == x

    # matrix containing non-integer rationals, small numerator and denominator
    x = Nemo.QQ[fmpq(1, 2) 2; 3 4]
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test convert(fmpq_mat, val) == x
    @test fmpq_mat(val) == x

    # matrix containing non-integer rationals, large numerator and denominator
    x = Nemo.QQ[fmpq(fmpz(2)^65, fmpz(3)^40) 2; 3 4]
    val = GAP.evalstr( "[ [ 2^65/3^40, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test convert(fmpq_mat, val) == x
    @test fmpq_mat(val) == x

    # matrix containing non-rationals
    val = GAP.evalstr( "[ [ E(4), 2 ], [ 3, 4 ] ]" )
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpq_mat, val)
end

@testset "matrices over a cyclotomic field" begin
    g = small_group(64, 140)
    reps = GAP.Globals.IrreducibleRepresentations(g.X)
    gmats = GAP.Globals.GeneratorsOfGroup(GAP.Globals.Image(reps[end]))
    omats, F, z = Oscar.matrices_over_cyclotomic_field(gmats)

    @test all(isone, [x^4 for x in omats])
end
