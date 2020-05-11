using Test, StraightLinePrograms, AbstractAlgebra

using StraightLinePrograms: Const, ExpPoly, Gen, MinusPoly, PlusPoly, RecPoly,
    TimesPoly, UniMinusPoly, pushconst!, pushop!,
    tmpmark, Line, Arg

const SL = StraightLinePrograms

replstr(c) = sprint((io, x) -> show(io, "text/plain", x), c)


@testset "LazyPolyRing" begin
    F = LazyPolyRing(ZZ)
    @test F isa LazyPolyRing{elem_type(ZZ)}
    @test F isa MPolyRing{elem_type(ZZ)}
    @test base_ring(F) == ZZ
end

@testset "RecPoly" begin
    # Const
    c = Const(1)
    @test c isa Const{Int}
    @test c isa RecPoly{Int}
    @test string(c) == "1"

    # Gen
    g = Gen{Int}(:x)
    @test g isa Gen{Int}
    @test g isa RecPoly{Int}
    @test string(g) == "x"

    # PlusPoly
    p = PlusPoly(c, g)
    @test p isa PlusPoly{Int} <: RecPoly{Int}
    @test p.xs[1] == c && p.xs[2] == g
    @test string(p) == "(1 + x)"

    # MinusPoly
    m = MinusPoly(p, g)
    @test m isa MinusPoly{Int} <: RecPoly{Int}
    @test string(m) == "((1 + x) - x)"

    # UniMinus
    u = UniMinusPoly(p)
    @test u isa UniMinusPoly{Int} <: RecPoly{Int}
    @test string(u) == "(-(1 + x))"

    # TimesPoly
    t = TimesPoly(g, p)
    @test t isa TimesPoly{Int} <: RecPoly{Int}
    @test string(t) == "(x(1 + x))"

    # ExpPoly
    e = ExpPoly(p, 3)
    @test e isa ExpPoly{Int} <: RecPoly{Int}
    @test string(e) == "(1 + x)^3"

    # +
    p1 =  e + t
    @test p1 isa PlusPoly{Int}
    @test p1.xs[1] === e
    @test p1.xs[2] === t
    p2 = p + e
    @test p2 isa PlusPoly{Int}
    @test p2.xs[1] === p.xs[1]
    @test p2.xs[2] === p.xs[2]
    @test p2.xs[3] === e
    p3 = e + p
    @test p3 isa PlusPoly{Int}
    @test p3.xs[1] === e
    @test p3.xs[2] === p.xs[1]
    @test p3.xs[3] === p.xs[2]
    p4 = p + p
    @test p4 isa PlusPoly{Int}
    @test p4.xs[1] === p.xs[1]
    @test p4.xs[2] === p.xs[2]
    @test p4.xs[3] === p.xs[1]
    @test p4.xs[4] === p.xs[2]

    # -
    m1 = e - t
    @test m1 isa MinusPoly
    @test m1.p === e
    @test m1.q === t
    m2 = -e
    @test m2 isa UniMinusPoly
    @test m2.p === e

    # *
    t1 =  e * p
    @test t1 isa TimesPoly{Int}
    @test t1.xs[1] === e
    @test t1.xs[2] === p
    t2 = t * e
    @test t2 isa TimesPoly{Int}
    @test t2.xs[1] === t.xs[1]
    @test t2.xs[2] === t.xs[2]
    @test t2.xs[3] === e
    t3 = e * t
    @test t3 isa TimesPoly{Int}
    @test t3.xs[1] === e
    @test t3.xs[2] === t.xs[1]
    @test t3.xs[3] === t.xs[2]
    t4 = t * t
    @test t4 isa TimesPoly{Int}
    @test t4.xs[1] === t.xs[1]
    @test t4.xs[2] === t.xs[2]
    @test t4.xs[3] === t.xs[1]
    @test t4.xs[4] === t.xs[2]

    # adhoc *
    am1 = 3 * p
    @test am1 isa RecPoly{Int}
    am2 = big(3) * p
    @test am2 isa RecPoly{Int}
    am3 = p * 3
    @test am3 isa RecPoly{Int}
    am4 = p * big(3)
    @test am4 isa RecPoly{Int}

    # ^
    e1 = p^3
    @test e1 isa ExpPoly
    @test e1.p === p
    @test e1.e == 3
end

@testset "LazyPoly" begin
    F = LazyPolyRing(zz)
    r = Const(1) + Gen{Int}(:x)
    p = LazyPoly(F, r)
    @test parent(p) === F
    @test string(p) == "(1 + x)"
    for x in (gen(F, :x), F(:x))
        @test x isa LazyPoly{Int}
        @test x.p isa Gen{Int}
        @test x.p.g == :x
    end
    c1 = F(2)
    @test c1 isa LazyPoly{Int}
    @test c1.p isa Const{Int}
    @test c1.p.c === 2

    @test (p+c1).p isa PlusPoly
    @test (p-c1).p isa MinusPoly
    @test (-p).p isa UniMinusPoly
    @test (p*c1).p isa TimesPoly
    @test (p^3).p isa ExpPoly
    @test_throws ArgumentError LazyPolyRing(ZZ)(big(1)) + c1
end

@testset "SLPolyRing" begin
    S = SLPolyRing(zz, [:x, :y])
    @test S isa SLPolyRing{Int}
    @test base_ring(S) == zz
    @test symbols(S) == [:x, :y]

    R, (x1, y1) = PolynomialRing(zz, ["x", "y"])

    x0, y0 = gens(S)
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
end

@testset "SLPoly" begin
    S = SLPolyRing(zz, [:x, :y])
    p = SLPoly(S, Int[], UInt64[])
    @test p isa SLPoly{Int,typeof(S)} <: MPolyElem{Int}
    @test parent(p) === S
    p = SLPoly(S)
    @test p isa SLPoly{Int,typeof(S)} <: MPolyElem{Int}
    @test parent(p) === S

    # copy
    q = SLPoly(S, Int[], UInt64[])
    # TODO: do smthg more interesting with q
    push!(q.cs, 3)
    push!(q.lines, Line(0))
    copy!(p, q)
    p2 = copy(q)
    for p1 in (p, p2)
        @test p1.cs == q.cs && p1.cs !== q.cs
        @test p1.lines == q.lines && p1.lines !== q.lines
    end
    S2 = SLPolyRing(zz, [:z, :t])
    @test_throws ArgumentError copy!(SLPoly(S2, Int[], UInt64[]), p)

    # building
    p = SLPoly(S)
    l1 = pushconst!(p, 1)
    @test p.cs == [1]
    @test l1 === Arg(1)
    l2 = pushconst!(p, 3)
    @test p.cs == [1, 3]
    @test l2 === Arg(2)

    l3 = pushop!(p, SL.plus, l1, l2)
    @test l3 == Arg(tmpmark | UInt64(1))
    @test p.lines[1].x == 0x0300000010000002
    l4 = pushop!(p, SL.times, l3, SL.input(2))
    @test l4 == Arg(tmpmark | UInt64(2))
    @test p.lines[2].x == 0x0540000018000002
    lines = copy(p.lines)
    @test p === SL.pushfinalize!(p, l4)
    @test p.lines == [Line(0x0300000010000002), Line(0x0500000038000002)]
    SL.pushinit!(p)
    @test lines == p.lines
    SL.pushfinalize!(p, l4)
    @test p.lines == [Line(0x0300000010000002), Line(0x0500000038000002)]
    # p == (1+3)*y
    @test SL.evaluate!(Int[], p, [1, 2]) == 8
    @test SL.evaluate!(Int[], p, [0, 3]) == 12
    l5 = SL.pushinit!(p)
    l6 = SL.pushop!(p, SL.times, SL.input(1), SL.input(2)) # xy
    l7 = SL.pushop!(p, SL.exponentiate, l5, Arg(2)) # (4y)^2
    l8 = SL.pushop!(p, SL.minus, l6, l7) # xy - 16y^2
    SL.pushfinalize!(p, l8)
    @test string(p) == "((xy) - ((1 + 3)y)^2)"
    @test SL.evaluate!(Int[], p, [2, 3]) == -138
    @test SL.evaluate!(Int[], p, [-2, -1]) == -14

    @test SL.evaluate(p, [2, 3]) == -138
    @test SL.evaluate(p, [-2, -1]) == -14

    # compile!
    pf = SL.compile!(p)
    @test pf([2, 3]) == -138
    @test pf([-2, -1]) == -14
    res = Int[]
    for xy in eachcol(rand(-99:99, 2, 100))
        v = Vector(xy) # TODO: don't require this
        @test pf(v) == SL.evaluate(p, v) == SL.evaluate!(res, p, v)
    end

    # conversion -> MPoly
    R, (x1, y1) = PolynomialRing(zz, ["x", "y"])
    q = convert(R, p)
    @test q isa Generic.MPoly
    @test parent(q) === R
    @test q == x1*y1 - 16*y1^2
    R2, (x2, y2) = PolynomialRing(zz, ["y", "x"])
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

    # construction from LazyPoly
    L = LazyPolyRing(zz)
    x, y = L(:x), L(:y)
    q = S(L(1))
    @test string(q) == "1"
    @test convert(R, q) == R(1)
    @test convert(R, S(x*y^2-x)) == x1*y1^2-x1
    @test convert(R, S(-(x+2*y)^3-4)) == -(x1+2*y1)^3-4

    # corner cases
    @test_throws BoundsError convert(R, SLPoly(S)) # error can change
    @test convert(R, S(L(1))) == R(1)
    @test convert(R, S(L(:x))) == x1

    # mutating ops
    # currently: p == xy - 16y^2
    # SL.addeq!(p, S(x)) # TODO: this bugs
    @test p === SL.addeq!(p, S(x*y))
    @test convert(R, p) == 2*x1*y1-16y1^2
    @test p === SL.subeq!(p, S(-16y^2))
    @test convert(R, p) == 2*x1*y1
    @test p === SL.subeq!(p)
    @test convert(R, p) == -2*x1*y1
    # @test p === SL.muleq!(p, p) # TODO: this bugs
    @test p === SL.muleq!(p, S(-2*x*y))
    @test convert(R, p) == 4*(x1*y1)^2
    @test p === SL.expeq!(p, 3)
    @test convert(R, p) == 64*(x1*y1)^6

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
end

@testset "SL internals" begin
    @test SL.showop == Dict(SL.uniplus      => '+',
                            SL.plus         => '+',
                            SL.uniminus     => '-',
                            SL.minus        => '-',
                            SL.times        => '*',
                            SL.divide       => '/',
                            SL.exponentiate => '^')
    @test length(SL.showop) == 7 # tests all keys are distinct
    for op in keys(SL.showop)
        @test SL.isuniplus(op) == (op == SL.uniplus)
        @test SL.istimes(op) == (op == SL.times)
        # ...
        @test (op.x & 0x8000000000000000 != 0) ==
            SL.isquasiunary(op) ==
            (op ∈ (SL.uniplus, SL.uniminus, SL.exponentiate))
        @test SL.isunary(op) == (op ∈ (SL.uniplus, SL.uniminus))
    end

    # pack & unpack
    ops = SL.Op.(rand(UInt64(0):UInt64(0xff), 100) .<< 62)
    is = rand(UInt64(0):SL.argmask, 100)
    js = rand(UInt64(0):SL.argmask, 100)
    @test SL.unpack.(SL.pack.(ops, is, js)) == tuple.(ops, Arg.(is), Arg.(js))

    for x = rand(Int64(0):Int(SL.tmpmark-1), 100)
        if SL.isinput(Arg(x))
            @test SL.input(x) == Arg(x)
        else
            @test SL.input(x).x ⊻ SL.inputmark == x
        end
    end
end
