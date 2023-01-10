@testset "fmpz" begin
    # small (GAP) integer
    x = fmpz(17)
    val = 17
    @test GAP.gap_to_julia(fmpz, val) == x
    @test fmpz(val) == x
    @test ZZ(val) == x

    # large positive GAP integer
    x = fmpz(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.gap_to_julia(fmpz, val) == x
    @test fmpz(val) == x
    @test ZZ(val) == x

    # large negative GAP integer
    x = -fmpz(2)^65
    val = GAP.evalstr("-2^65")
    @test GAP.gap_to_julia(fmpz, val) == x
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
    @test fmpq(val) == x
    @test QQ(val) == x

    # large positive GAP integer
    x = fmpq(2)^65
    val = GAP.evalstr("2^65")
    @test GAP.gap_to_julia(fmpq, val) == x
    @test fmpq(val) == x
    @test QQ(val) == x

    # large negative GAP integer
    x = -fmpq(2)^65
    val = GAP.evalstr("-2^65")
    @test GAP.gap_to_julia(fmpq, val) == x
    @test fmpq(val) == x
    @test QQ(val) == x

    # s "proper" rationals with large and small numerators and denominators
    @testset "fmpq $a / $b" for a in [2, -2, fmpz(2^65), -fmpz(2^65)], b in [3, -3, fmpz(3^40), -fmpz(3^50)]
        x = fmpq(a, b)
        val = GAP.evalstr("$a/$b")
        @test GAP.gap_to_julia(fmpq, val) == x
        @test fmpq(val) == x
        @test QQ(val) == x
    end

    # non-rational
    val = GAP.evalstr("()")
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpq, val)
end

@testset "fmpz_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.ZZ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpz_mat, val) == x
    @test fmpz_mat(val) == x
    @test matrix(ZZ, val) == x

    # matrix containing small and large integers
    x = Nemo.ZZ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpz_mat, val) == x
    @test fmpz_mat(val) == x
    @test matrix(ZZ, val) == x

    # matrix containing non-integers
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpz_mat, val)
end

@testset "fmpq_mat" begin
    # matrix of small (GAP) integers
    x = Nemo.QQ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test fmpq_mat(val) == x
    @test matrix(QQ, val) == x

    # matrix containing small and large integers
    x = Nemo.QQ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test fmpq_mat(val) == x
    @test matrix(QQ, val) == x

    # matrix containing non-integer rationals, small numerator and denominator
    x = Nemo.QQ[fmpq(1, 2) 2; 3 4]
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test fmpq_mat(val) == x
    @test matrix(QQ, val) == x

    # matrix containing non-integer rationals, large numerator and denominator
    x = Nemo.QQ[fmpq(fmpz(2)^65, fmpz(3)^40) 2; 3 4]
    val = GAP.evalstr( "[ [ 2^65/3^40, 2 ], [ 3, 4 ] ]" )
    @test GAP.gap_to_julia(fmpq_mat, val) == x
    @test fmpq_mat(val) == x
    @test matrix(QQ, val) == x

    # matrix containing non-rationals
    val = GAP.evalstr( "[ [ E(4), 2 ], [ 3, 4 ] ]" )
    @test_throws GAP.ConversionError GAP.gap_to_julia(fmpq_mat, val)
end

@testset "finite field elements" begin
  @testset "with characteristic $p" for p in [ 5, fmpz(5), 65537, fmpz(65537) ]
    @testset "with finite field $F" for F in [ GF(p), GF(p,1), GF(p,2) ]
        x = F(2)
        q = order(F)

        # GAP large integers
        val = GAP.evalstr( "2+$p*2^65" )
        @test F(val) == x

        # finite field elements from prime field
        z_gap = GAP.evalstr( "Z($p)" )
        @test F(z_gap)^2 == F(z_gap^2)

        # finite field elements from full field
        # FIXME: this is not yet implemented in general
        #z_gap = GAP.evalstr( "Z($q)" )
        #@test F(z_gap)^2 == F(z_gap^2)
    end
  end
end

@testset "finite field matrix" begin
  @testset "with characteristic $p" for p in [ 5, fmpz(5), 65537, fmpz(65537) ]
    @testset "with finite field $F" for F in [ GF(p), GF(p,1), GF(p,2) ]
        x = F[1 2; 3 4]

        # matrix of small (GAP) integers
        val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
        @test matrix(F, val) == x

        # matrix containing small and large integers
        val = GAP.evalstr( "[ [ 1, 2+$p*2^65 ], [ 3, 4 ] ]" )
        @test matrix(F, val) == x

        # matrix of finite field elements
        val = GAP.evalstr( "Z($p)^0 * [ [ 1, 2 ], [ 3, 4 ] ]" )
        @test matrix(F, val) == x

        # possible compressed matrix of finite field elements
        GAP.Globals.ConvertToMatrixRep(val)
        @test matrix(F, val) == x
    end
  end
end

@testset "single cyclotomics" begin
    # to cyclotomic fields
    F, z = CyclotomicField(1)
    @test F(GAP.evalstr("2^64")) == F(2)^64

    F, z = CyclotomicField(5)
    @test F(GAP.Globals.Sqrt(5)) == -2*z^3 - 2*z^2 - 1

    F, z = CyclotomicField(15)
    @test F(GAP.Globals.E(5)) == z^3
    @test F(GAP.Globals.E(3)) == z^5

    # to `QQAbElem`
    x = QQAbElem(GAP.evalstr("2^64"))
    @test x == fmpz(2)^64

    F, z = abelian_closure(QQ)
    val = GAP.evalstr("EB(5)")
    x = QQAbElem(val)
    @test x == z(5) + z(5)^4
    @test F(val) == x
    @test GAP.gap_to_julia(QQAbElem, val) == x

    # not supported conversions
    F, z = quadratic_field(5)
    @test_throws ArgumentError F(GAP.Globals.Sqrt(5))

    F, z = CyclotomicField(5)
    @test_throws ArgumentError F(GAP.Globals.Sqrt(7))

    F, z = CyclotomicField(5)
    x = GAP.Globals.Indeterminate(GAP.Globals.Rationals)
    pol = x^2 - 5
    gapF = GAP.Globals.AlgebraicExtension(GAP.Globals.Rationals, pol)
    a = GAP.Globals.PrimitiveElement(gapF)
    @test_throws ArgumentError F(a)

    @test_throws GAP.ConversionError QQAbElem(GAP.evalstr("[ E(3) ]"))
end

@testset "matrices over a cyclotomic field" begin
    g = small_group(64, 140)
    reps = GAP.Globals.IrreducibleRepresentations(g.X)
    gmats = GAP.Globals.GeneratorsOfGroup(GAP.Globals.Image(reps[end]))
    omats, F, z = Oscar.matrices_over_cyclotomic_field(gmats)

    @test all(isone, [x^4 for x in omats])
    @test all(i -> omats[i] == matrix(F, gmats[i]), 1:length(gmats))

    F, _ = CyclotomicField(4)
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(F, GAP.Globals.E(4))
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(F, GAP.evalstr("[]"))
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(F, GAP.evalstr("[ [ [ Z(2) ] ] ]"))
    @test Oscar.matrices_over_cyclotomic_field(F, gmats) == (omats, F, z)

    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(GAP.Globals.E(4))
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(GAP.evalstr("[]"))
    @test_throws MethodError Oscar.matrices_over_cyclotomic_field(GAP.evalstr("[ [ [ Z(2) ] ] ]"))

    F, z = quadratic_field(5)
    @test_throws ArgumentError matrix(F, GAP.evalstr("[ [ Sqrt(5) ] ]"))

    F, z = CyclotomicField(5)
    x = GAP.Globals.Indeterminate(GAP.Globals.Rationals)
    pol = x^2 - 5
    gapF = GAP.Globals.AlgebraicExtension(GAP.Globals.Rationals, pol)
    a = GAP.Globals.PrimitiveElement(gapF)
    @test_throws ArgumentError F(GAP.Globals.IdentityMat(2)*a)

    F, z = CyclotomicField(15)
    @test matrix(F, GAP.evalstr("[ [ Sqrt(5) ] ]"))[1,1] == -2*z^7 + 2*z^3 - 2*z^2 + 1

    F, z = CyclotomicField(5)
    @test matrix(F, GAP.evalstr("[ [ E(5) ] ]"))[1,1] == z

    F, z = CyclotomicField(7)
    @test_throws ErrorException matrix(F, GAP.evalstr("[ [ E(5) ] ]"))
end

@testset "matrices over a field" begin
    # internal finite field elements in GAP
    mats = GAP.evalstr("[ [ [ Z(2) ] ], [ [ Z(4) ] ] ]")
    omats, F, z = Oscar.matrices_over_cyclotomic_field(mats)
    @test order(F) == 4
    @test omats[1][1,1] == one(F)
    @test omats[2][1,1] == gen(F)

    # algebraic extension elements in GAP
    oF = GF(2)
    gF = Oscar.iso_oscar_gap(oF)
    R, x = PolynomialRing(oF)
    gpol = image(Oscar.iso_oscar_gap(R), x^2+x+1)
    F = GAP.Globals.AlgebraicExtension(codomain(gF), gpol)
    mats = GAP.GapObj([[[GAP.Globals.One(F)]], [GAP.Globals.GeneratorsOfField(F)]], recursive = true)
    omats, F, z = Oscar.matrices_over_cyclotomic_field(mats)
    @test order(F) == 4
    @test omats[1][1,1] == one(F)
    @test omats[2][1,1] == gen(F)

    # not a collection in GAP
    mats = GAP.evalstr("[ [ [ Z(2) ] ], [ [ Z(3) ] ] ]")
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(mats)
end
