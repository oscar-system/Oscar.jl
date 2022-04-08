using Oscar: SLPolynomialRing

replstr(c) = sprint((io, x) -> show(io, "text/plain", x), c)

@testset "LazyPolyRing" begin
    F = LazyPolyRing(AbstractAlgebra.ZZ)
    @test F isa LazyPolyRing{elem_type(AbstractAlgebra.ZZ)}
    @test F isa MPolyRing{elem_type(AbstractAlgebra.ZZ)}
    @test base_ring(F) == AbstractAlgebra.ZZ
end


@testset "LazyPoly" begin
    F = LazyPolyRing(AbstractAlgebra.zz)
    r = SLP.Const(1) + SLP.Gen(:x)
    p = LazyPoly(F, r)
    @test parent(p) === F
    @test string(p) == "(1 + x)"
    for x in (gen(F, :x), F(:x))
        @test x isa LazyPoly{Int}
        @test x.p isa SLP.Gen
        @test x.p.g == :x
    end
    c1 = F(2)
    @test c1 isa LazyPoly{Int}
    @test c1.p isa SLP.Const{Int}
    @test c1.p.c === 2

    @test (p+c1).p isa SLP.Plus
    @test (p-c1).p isa SLP.Minus
    @test (-p).p isa SLP.UniMinus
    @test (p*c1).p isa SLP.Times
    #@test (p^3).p isa SLP.Exp
    @test_throws ArgumentError LazyPolyRing(AbstractAlgebra.ZZ)(big(1)) + c1
end

@testset "SLPolyRing" begin
    S = SLPolyRing(AbstractAlgebra.zz, [:x, :y])
    @test S isa SLPolyRing{Int}
    @test base_ring(S) == AbstractAlgebra.zz
    @test symbols(S) == [:x, :y]

    for S2 in (SLPolyRing(AbstractAlgebra.zz, ["x", "y"]),
               SLPolyRing(AbstractAlgebra.zz, ['x', 'y']))
        @test S2 isa SLPolyRing{Int}
        @test base_ring(S2) == AbstractAlgebra.zz
        @test symbols(S2) == [:x, :y]
    end
    for Sxy in (SLPolynomialRing(AbstractAlgebra.zz, ["x", "y"]),
                SLPolynomialRing(AbstractAlgebra.zz, ['x', 'y']))
        S3, (x, y) = Sxy
        @test S3 isa SLPolyRing{Int}
        @test base_ring(S3) == AbstractAlgebra.zz
        @test symbols(S3) == [:x, :y]
        @test string(x) == "x"
        @test string(y) == "y"
        @test parent(x) == S3
        @test parent(y) == S3
    end

    S4 = SLPolyRing(AbstractAlgebra.zz, 3)
    @test S4 isa SLPolyRing{Int}
    @test ngens(S4) == 3
    X = gens(S4)
    @test length(X) == 3
    @test string.(X) == ["x1", "x2", "x3"]

    S4, X = SLPolynomialRing(AbstractAlgebra.zz, 0x2)
    @test S4 isa SLPolyRing{Int}
    @test ngens(S4) == 2
    X = gens(S4)
    @test length(X) == 2
    @test string.(X) == ["x1", "x2"]

    S5 = SLPolyRing(AbstractAlgebra.zz, :x => 1:3, :y => [2, 4])
    @test S5 isa SLPolyRing{Int}
    XS = gens(S5)
    @test string.(XS) == ["x1", "x2", "x3", "y2", "y4"]

    S5, X, Y = SLPolynomialRing(AbstractAlgebra.zz, :x => 1:3, "y" => [2, 4])
    @test S5 isa SLPolyRing{Int}
    XS = gens(S5)
    @test string.(XS) == ["x1", "x2", "x3", "y2", "y4"]
    @test X == XS[1:3]
    @test Y == XS[4:5]

    s1 = one(S)
    @test s1 == 1
    @test s1 isa SLPoly{Int}
    s0 = zero(S)
    @test s0 == 0
    @test s0 isa SLPoly{Int}

    R, (x1, y1) = PolynomialRing(AbstractAlgebra.zz, ["x", "y"])

    x0, y0 = gens(S)
    @test ngens(S) == 2
    @test nvars(S) == 2
    @test string(x0) == "x"
    @test string(y0) == "y"
    @test replstr(x0) == "x"
    @test replstr(y0) == "y"
    @test string(S(2)) == "2"
    @test replstr(S(2)) == "2"

    for x in (gen(S, 1), x0)
        @test string(x) == "x"
        @test x isa SLPoly{Int}
        @test convert(R, x) == x1
    end
    for y in (gen(S, 2), y0)
        @test string(y) == "y"
        @test y isa SLPoly{Int}
        @test convert(R, y) == y1
    end

    for t = (2, big(2), 0x2)
        @test S(t) isa SLPoly{Int,typeof(S)}
    end
end

@testset "SLPoly" begin
    S = SLPolyRing(AbstractAlgebra.zz, [:x, :y])
    p = SLPoly(S, SLP.SLProgram{Int}())
    @test p isa SLPoly{Int,typeof(S)}
    @test parent(p) === S
    @test parent_type(p) == parent_type(typeof(p)) == typeof(S)
    @test elem_type(S) == elem_type(typeof(S)) == typeof(p)

    p = SLPoly(S)
    @test p isa SLPoly{Int,typeof(S)}
    @test parent(p) === S

    @test S(p) === p

    @test zero(p) == zero(S)
    @test zero(p) isa SLPoly{Int}
    @test one(p) == one(S)
    @test one(p) isa SLPoly{Int}

    @test !isone(p)
    @test !iszero(p)

    # copy
    q = SLPoly(S)
    # TODO: do smthg more interesting with q
    push!(SLP.constants(q), 3)
    push!(q.slprogram.lines, SLP.Line(0))
    copy!(p, q)
    p2 = copy(q)
    for p1 in (p, p2)
        @test SLP.constants(p1) == SLP.constants(q) && SLP.constants(p1) !== SLP.constants(q)
        @test SLP.lines(p1) == SLP.lines(q) && SLP.lines(p1) !== SLP.lines(q)
    end
    S2 = SLPolyRing(AbstractAlgebra.zz, [:z, :t])
    @test_throws ArgumentError copy!(SLPoly(S2, SLP.SLProgram{Int}()), p)
    @test_throws ArgumentError S2(p) # wrong parent

    # building
    p = SLPoly(S)
    l1 = SLP.pushconst!(p.slprogram, 1)
    @test SLP.constants(p) == [1]
    @test l1 === SLP.asconstant(1)

    # currently not supported anymore
    # l2 = pushconst!(p, 3)
    # @test constants(p) == [1, 3]
    # @test l2 === SLP.asconstant(2)

    l3 = SLP.pushop!(p, SLP.plus, l1, SLP.input(1))
    @test l3 == SLP.Arg(UInt64(1))
    @test SLP.lines(p)[1].x == 0x0340000018000001
    l4 = SLP.pushop!(p, SLP.times, l3, SLP.input(2))
    @test l4 == SLP.Arg(UInt64(2))
    @test SLP.lines(p)[2].x ==0x0500000018000002
    pl = copy(SLP.lines(p))
    @test p === SLP.pushfinalize!(p, l4)

    @test SLP.lines(p) == [SLP.Line(0x0340000018000001), SLP.Line(0x0500000018000002)]
    SLP.pushinit!(p)
    @test pl == SLP.lines(p)
    SLP.pushfinalize!(p, l4)
    @test SLP.lines(p) == [SLP.Line(0x0340000018000001), SLP.Line(0x0500000018000002)]
    # p == (1+x)*y
    @test SLP.evaluate!(Int[], p, [1, 2]) == 4
    @test SLP.evaluate!(Int[], p, [0, 3]) == 3
    l5 = SLP.pushinit!(p)
    l6 = SLP.pushop!(p, SLP.times, SLP.input(1), SLP.input(2)) # xy
    l7 = SLP.pushop!(p, SLP.exponentiate, l5, SLP.Arg(2)) # ((1+x)y)^2
    l8 = SLP.pushop!(p, SLP.minus, l6, l7) # xy - ((1+x)y)^2
    SLP.pushfinalize!(p, l8)
    @test string(p) == "x*y - ((1 + x)*y)^2"
    @test SLP.evaluate!(Int[], p, [2, 3]) == -75
    @test SLP.evaluate!(Int[], p, [-2, -1]) == 1

    @test evaluate(p, [2, 3]) == -75
    @test evaluate(p, [-2, -1]) == 1

    # nsteps
    @test SLP.nsteps(p) == 5

    # compile!
    pf = SLP.compile!(p)
    @test pf([2, 3]) == -75
    @test pf([-2, -1]) == 1
    res = Int[]
    for xy in eachcol(rand(-99:99, 2, 100))
        v = Vector(xy) # TODO: don't require this
        @test pf(v) == evaluate(p, v) == SLP.evaluate!(res, p, v)
    end

    # conversion -> MPoly
    R, (x1, y1) = PolynomialRing(AbstractAlgebra.zz, ["x", "y"])
    q = convert(R, p)
    @test q isa Generic.MPoly
    @test parent(q) === R
    @test q == -x1^2*y1^2-2*x1*y1^2+x1*y1-y1^2
    R2, (x2, y2) = PolynomialRing(AbstractAlgebra.zz, ["y", "x"])
    @test_throws ArgumentError convert(R2, p)

    @test convert(R, S()) == R()
    @test convert(R, S(3)) == R(3)
    @test convert(R, S(-4)) == R(-4)

    # conversion MPoly -> SLPoly
    for _=1:100
        r = rand(R, 1:20, 0:13, -19:19)
        @test convert(R, convert(S, r)) == r
        @test convert(R, convert(S, r; limit_exp=true)) == r
    end
    r = R()
    @test convert(R, convert(S, r)) == r

    # construction from LazyPoly
    L = LazyPolyRing(AbstractAlgebra.zz)
    x, y = L(:x), L(:y)
    q = S(L(1))
    @test string(q) == "1"
    @test convert(R, q) == R(1)
    @test convert(R, S(x*y^2-x)) == x1*y1^2-x1
    @test convert(R, S(-(x+2*y)^3-4)) == -(x1+2*y1)^3-4

    # corner cases
    @test_throws ArgumentError convert(R, SLPoly(S)) # error can change
    @test convert(R, S(L(1))) == R(1)
    @test convert(R, S(L(:x))) == x1

    # mutating ops
    X, Y = gens(S)
    p = S(x*y-16*y^2)
    # SLP.addeq!(p, S(x)) # TODO: this bugs
    @test p === addeq!(p, S(x*y))
    @test convert(R, p) == 2*x1*y1-16y1^2
    @test p == X*Y-16Y^2+X*Y

    @test p === SLP.subeq!(p, S(-16y^2))
    @test convert(R, p) == 2*x1*y1
    @test p == X*Y-16Y^2+X*Y- (-16*Y^2)

    @test p === SLP.subeq!(p)
    @test convert(R, p) == -2*x1*y1
    # @test p === SLP.muleq!(p, p) # TODO: this bugs
    @test p === SLP.muleq!(p, S(-2*x*y))
    @test convert(R, p) == 4*(x1*y1)^2
    @test p === SLP.expeq!(p, 3)
    @test convert(R, p) == 64*(x1*y1)^6

    # permutegens! and ^
    _, (X1, X2, X3) = SLPolynomialRing(AbstractAlgebra.zz, [:x1, :x2, :x3])
    p = X1*X2^2+X3^3
    p0 = copy(p)
    perm = [3, 1, 2]
    q = Oscar.permutegens!(p, perm)
    @test q === p
    @test q == X3*X1^2+X2^3
    @test q == p0^Perm(perm)
    @test p0 == X1*X2^2+X3^3 # not mutated

    # binary/unary ops
    p = S(x*y - 16y^2)
    p = p + S(x*y)
    @test p isa SLPoly{Int}
    @test convert(R, p) == 2*x1*y1-16y1^2
    p = p - S(-16y^2)
    @test p isa SLPoly{Int}
    @test convert(R, p) == 2*x1*y1
    p = -p
    @test p isa SLPoly{Int}
    @test convert(R, p) == -2*x1*y1
    p = p * S(-2*x*y)
    @test p isa SLPoly{Int}
    @test convert(R, p) == 4*(x1*y1)^2
    p = p^3
    @test p isa SLPoly{Int}
    @test convert(R, p) == 64*(x1*y1)^6

    # adhoc ops
    p = S(x*y - 16y^2)
    q = convert(R, p)
    @test convert(R, 2p) == 2q
    @test convert(R, p*3) == q*3
    @test convert(R, 2+p) == 2+q
    @test convert(R, p+2) == q+2
    @test convert(R, 2-p) == 2-q
    @test convert(R, p-2) == q-2

    R = ResidueRing(AbstractAlgebra.ZZ, 3)
    a = R(2)
    S = SLPolyRing(R, [:x, :y])
    x, y = gens(S)
    @test parent(a*x) == S
    @test parent(x*a) == S
    @test parent(a+x) == S
    @test parent(x+a) == S
    @test parent(a-x) == S
    @test parent(x-a) == S

    # 3-args evaluate
    S = SLPolyRing(AbstractAlgebra.zz, [:x, :y])
    x, y = gens(S)
    p = 2*x^3+y^2+3
    @test evaluate(p, [2, 3]) == 28
    @test SLP.evaluate!(Int[], p, [2, 3]) == 28
    @test SLP.evaluate!(Int[], p, [2, 3], x -> x) == 28
    @test SLP.evaluate!(Int[], p, [2, 3], x -> -x) == -10

    # trivial rings
    S = SLPolyRing(AbstractAlgebra.zz, Symbol[])
    gs = gens(S)
    @test evaluate(S(1), gs) == S(1)
    @test SLP.evaluate!(empty(gs), S(1), gs) == S(1)

    # evaluate MPoly at SLPolyRing generators
    R, (x, y) = PolynomialRing(AbstractAlgebra.zz, ["x", "y"])
    S = SLPolyRing(AbstractAlgebra.zz, [:x, :y])
    X, Y = gens(S)
    p = evaluate(x+y, [X, Y])
    # this is bad to hardcode exactly how evaluation of `x+y` happens,
    # we just want to test that this works and looks correct
    @test p ==  0 + 1*(1*X^1) + 1*(1*Y^1)

    @testset "SLPolyRing is a proper Ring" begin
        S, (x1, x2) = Oscar.SLPolynomialRing(QQ, 2)
        St, t = PolynomialRing(S, "t")

        @testset "show" begin
            @test string(x1) == "x1"
            @test string(-x1*QQ(2, 3)*x2^3) == "-x1*2//3*x2^3"
            @test string(x1-x2+x1) == "x1 - x2 + x1"
            @test string(2-x2) == "2 - x2"
            @test string(2t+3) == "2*t + 3"
        end

        q = x1+x2
        @test q === zero!(q)
        @test string(q) == "0" # can't test with iszero, which currently always return false

        @test q === mul!(q, x1, x2)
        @test string(q) == "x1*x2"
        @test q === add!(q, x1, x2)
        @test string(q) == "x1 + x2"

        r = prod(t-y for y = gens(S))
        @test string(r) == "t^2 + (-x2 - x1)*t + x1*x2"
    end

    # issue #250
    let (S, xs) = slpoly_ring(QQ, 2)
        f = S(1)
        x = evaluate(f, xs)
        @test parent(x) == parent(f)
        @test typeof(x) == typeof(f)
    end

    # issue #253
    let (S, (f,)) = slpoly_ring(ZZ, 1)
        Zx, x = PolynomialRing(ZZ)
        f1 = evaluate(f, [x(f)])
        f2 = evaluate(f1, [f1])
        @test parent(f1) == parent(f2)
    end
end
