@testset "ZZRingElem" begin
    # small (GAP) integer
    x = ZZRingElem(17)
    val = 17
    @test ZZRingElem(val) == x
    @test ZZ(val) == x

    # large positive GAP integer
    x = ZZRingElem(2)^65
    val = GAP.evalstr("2^65")
    @test ZZRingElem(val) == x
    @test ZZ(val) == x

    # large negative GAP integer
    x = -ZZRingElem(2)^65
    val = GAP.evalstr("-2^65")
    @test ZZRingElem(val) == x
    @test ZZ(val) == x

    # non-integer
    val = GAP.evalstr("1/2")
    @test_throws GAP.ConversionError ZZRingElem(val)
    @test_throws GAP.ConversionError ZZ(val)
end

@testset "QQFieldElem" begin
    # small (GAP) integer
    x = QQFieldElem(17)
    val = 17
    @test QQFieldElem(val) == x
    @test QQ(val) == x

    # large positive GAP integer
    x = QQFieldElem(2)^65
    val = GAP.evalstr("2^65")
    @test QQFieldElem(val) == x
    @test QQ(val) == x

    # large negative GAP integer
    x = -QQFieldElem(2)^65
    val = GAP.evalstr("-2^65")
    @test QQFieldElem(val) == x
    @test QQ(val) == x

    # s "proper" rationals with large and small numerators and denominators
    @testset "QQFieldElem $a / $b" for a in [2, -2, ZZRingElem(2^65), -ZZRingElem(2^65)], b in [3, -3, ZZRingElem(3^40), -ZZRingElem(3^50)]
        x = QQFieldElem(a, b)
        val = GAP.evalstr("$a/$b")
        @test QQFieldElem(val) == x
        @test QQ(val) == x
    end

    # non-rational
    val = GAP.evalstr("E(4)")
    @test_throws GAP.ConversionError QQFieldElem(val)
    @test_throws GAP.ConversionError QQ(val)
end

@testset "ZZMatrix" begin
    # matrix of small (GAP) integers
    x = Nemo.ZZ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test ZZMatrix(val) == x
    @test matrix(ZZ, val) == x

    # matrix containing small and large integers
    x = Nemo.ZZ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test ZZMatrix(val) == x
    @test matrix(ZZ, val) == x

    # matrix containing non-integers
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test_throws GAP.ConversionError ZZMatrix(val)
    @test_throws GAP.ConversionError matrix(ZZ, val)
end

@testset "QQMatrix" begin
    # matrix of small (GAP) integers
    x = Nemo.QQ[1 2; 3 4]
    val = GAP.evalstr( "[ [ 1, 2 ], [ 3, 4 ] ]" )
    @test QQMatrix(val) == x
    @test matrix(QQ, val) == x

    # matrix containing small and large integers
    x = Nemo.QQ[1 BigInt(2)^65; 3 4]
    val = GAP.evalstr( "[ [ 1, 2^65 ], [ 3, 4 ] ]" )
    @test QQMatrix(val) == x
    @test matrix(QQ, val) == x

    # matrix containing non-integer rationals, small numerator and denominator
    x = Nemo.QQ[QQFieldElem(1, 2) 2; 3 4]
    val = GAP.evalstr( "[ [ 1/2, 2 ], [ 3, 4 ] ]" )
    @test QQMatrix(val) == x
    @test matrix(QQ, val) == x

    # matrix containing non-integer rationals, large numerator and denominator
    x = Nemo.QQ[QQFieldElem(ZZRingElem(2)^65, ZZRingElem(3)^40) 2; 3 4]
    val = GAP.evalstr( "[ [ 2^65/3^40, 2 ], [ 3, 4 ] ]" )
    @test QQMatrix(val) == x
    @test matrix(QQ, val) == x

    # matrix containing non-rationals
    val = GAP.evalstr( "[ [ E(4), 2 ], [ 3, 4 ] ]" )
    @test_throws GAP.ConversionError QQMatrix(val)
    @test_throws GAP.ConversionError matrix(QQ, val)
end

@testset "finite field elements" begin
  @testset "with characteristic $p" for p in [ 5, ZZRingElem(5), 65537, ZZRingElem(65537) ]
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
  @testset "with characteristic $p" for p in [ 5, ZZRingElem(5), 65537, ZZRingElem(65537) ]
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
    F, z = cyclotomic_field(1)
    @test F(GAP.evalstr("2^64")) == F(2)^64

    F, z = cyclotomic_field(5)
    @test F(GAP.Globals.Sqrt(5)) == -2*z^3 - 2*z^2 - 1

    F, z = cyclotomic_field(15)
    @test F(GAP.Globals.E(5)) == z^3
    @test F(GAP.Globals.E(3)) == z^5

    # to `QQAbFieldElem`
    x = QQAbFieldElem(GAP.evalstr("2^64"))
    @test x == ZZRingElem(2)^64

    F, z = abelian_closure(QQ)
    val = GAP.evalstr("EB(5)")
    x = QQAbFieldElem(val)
    @test x == z(5) + z(5)^4
    @test F(val) == x
    @test QQAbFieldElem(val) == x

    # not supported conversions
    F, z = quadratic_field(5)
    @test_throws ArgumentError F(GAP.Globals.Sqrt(5))

    F, z = cyclotomic_field(5)
    @test_throws ArgumentError F(GAP.Globals.Sqrt(7))

    F, z = cyclotomic_field(5)
    x = GAP.Globals.Indeterminate(GAP.Globals.Rationals)
    pol = x^2 - 5
    gapF = GAP.Globals.AlgebraicExtension(GAP.Globals.Rationals, pol)
    a = GAP.Globals.PrimitiveElement(gapF)
    @test_throws ArgumentError F(a)

    @test_throws GAP.ConversionError QQAbFieldElem(GAP.evalstr("[ E(3) ]"))
end

@testset "matrices over a cyclotomic field" begin
    g = small_group(64, 140)
    reps = GAP.Globals.IrreducibleRepresentations(GapObj(g))
    gmats = GAP.Globals.GeneratorsOfGroup(GAP.Globals.Image(reps[end]))
    omats, F, z = Oscar.matrices_over_cyclotomic_field(gmats)

    @test all(isone, [x^4 for x in omats])
    @test all(i -> omats[i] == matrix(F, gmats[i]), 1:length(gmats))

    F, _ = cyclotomic_field(4)
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(F, GAP.Globals.E(4))
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(F, GAP.evalstr("[]"))
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(F, GAP.evalstr("[ [ [ Z(2) ] ] ]"))
    @test Oscar.matrices_over_cyclotomic_field(F, gmats) == (omats, F, z)

    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(GAP.Globals.E(4))
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(GAP.evalstr("[]"))

    F, z = quadratic_field(5)
    @test_throws ArgumentError matrix(F, GAP.evalstr("[ [ Sqrt(5) ] ]"))

    F, z = cyclotomic_field(5)
    x = GAP.Globals.Indeterminate(GAP.Globals.Rationals)
    pol = x^2 - 5
    gapF = GAP.Globals.AlgebraicExtension(GAP.Globals.Rationals, pol)
    a = GAP.Globals.PrimitiveElement(gapF)
    @test_throws ArgumentError F(GAP.Globals.IdentityMat(2)*a)

    F, z = cyclotomic_field(15)
    @test matrix(F, GAP.evalstr("[ [ Sqrt(5) ] ]"))[1,1] == -2*z^7 + 2*z^3 - 2*z^2 + 1

    F, z = cyclotomic_field(5)
    @test matrix(F, GAP.evalstr("[ [ E(5) ] ]"))[1,1] == z

    F, z = cyclotomic_field(7)
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
    R, x = polynomial_ring(oF)
    gpol = image(Oscar.iso_oscar_gap(R), x^2+x+1)
    F = GAP.Globals.AlgebraicExtension(codomain(gF), gpol)
    mats = GapObj([[[GAP.Globals.One(F)]], [GAP.Globals.GeneratorsOfField(F)]]; recursive = true)
    omats, F, z = Oscar.matrices_over_cyclotomic_field(mats)
    @test order(F) == 4
    @test omats[1][1,1] == one(F)
    @test omats[2][1,1] == gen(F)

    # not a collection in GAP
    mats = GAP.evalstr("[ [ [ Z(2) ] ], [ [ Z(3) ] ] ]")
    @test_throws ArgumentError Oscar.matrices_over_cyclotomic_field(mats)
end

@testset "straight line programs" begin
    l = GapObj([[[1, 2], 3], [[3, 2], 2], [1, 2, 2, 1]], recursive = true)
    gapslp = GAP.Globals.StraightLineProgram(l, 2)
    slp = straight_line_program(gapslp)
    @test evaluate(slp, [2, 3]) == 64
    @test GAP.Globals.ResultOfStraightLineProgram(gapslp, GapObj([2, 3])) == 64

    @test_throws ArgumentError straight_line_program(l)
end
