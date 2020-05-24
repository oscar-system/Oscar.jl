using Test, StraightLinePrograms, AbstractAlgebra

using StraightLinePrograms: Const, Exp, Gen, Minus, Plus, Lazy,
    Times, UniMinus, pushconst!, pushop!,
    Line, Arg, constants, lines, evaluate

const SL = StraightLinePrograms

replstr(c) = sprint((io, x) -> show(io, "text/plain", x), c)


@testset "LazyPolyRing" begin
    F = LazyPolyRing(ZZ)
    @test F isa LazyPolyRing{elem_type(ZZ)}
    @test F isa MPolyRing{elem_type(ZZ)}
    @test base_ring(F) == ZZ
end

@testset "Lazy" begin
    x, y, z = Gen.([:x, :y, :z])

    # Const
    c = Const(1)
    @test c isa Const{Int}
    @test c isa Lazy
    @test string(c) == "1"
    @test isempty(SL.gens(c))
    @test c == Const(1) == Const(0x1)
    @test c != Const(2)

    # Gen
    g = Gen(:x)
    @test g isa Gen
    @test g isa Lazy
    @test string(g) == "x"
    @test SL.gens(g) == [:x]
    @test g == x
    @test g != y

    # Plus
    p = Plus(c, g)
    @test p isa Plus <: Lazy
    @test p.xs[1] == c && p.xs[2] == g
    @test string(p) == "(1 + x)"
    @test SL.gens(p) == [:x]
    @test p == 1+x == 0x1+x
    @test p != 2+x && p != 1+y

    # Minus
    m = Minus(p, g)
    @test m isa Minus <: Lazy
    @test string(m) == "((1 + x) - x)"
    @test SL.gens(m) == [:x]
    @test m == (1+x)-x
    @test m != (1+x)+x && m != x-x && m != (1+x)-y

    # UniMinus
    u = UniMinus(p)
    @test u isa UniMinus <: Lazy
    @test string(u) == "(-(1 + x))"
    @test SL.gens(u) == [:x]
    @test u == -(1+x)
    @test u != (1+x) && u != -(1+y)

    # Times
    t = Times(g, p)
    @test t isa Times <: Lazy
    @test string(t) == "(x(1 + x))"
    @test SL.gens(t) == [:x]
    @test t == x*(1+x)
    @test t != (1+x)*x && t != y*(1+x) && t != x*(1+y)

    # Exp
    e = Exp(p, 3)
    @test e isa Exp <: Lazy
    @test string(e) == "(1 + x)^3"
    @test SL.gens(e) == [:x]
    @test e == (1+x)^3
    @test e != (1+x)^4 && e != (1+y)^3

    # +
    p1 =  e + t
    @test p1 isa Plus
    @test p1.xs[1] === e
    @test p1.xs[2] === t
    @test p1 == e+t
    p2 = p + e
    @test p2 isa Plus
    @test p2.xs[1] === p.xs[1]
    @test p2.xs[2] === p.xs[2]
    @test p2.xs[3] === e
    @test p2 == p+e
    p3 = e + p
    @test p3 isa Plus
    @test p3.xs[1] === e
    @test p3.xs[2] === p.xs[1]
    @test p3.xs[3] === p.xs[2]
    @test p3 == e+p
    p4 = p + p
    @test p4 isa Plus
    @test p4.xs[1] === p.xs[1]
    @test p4.xs[2] === p.xs[2]
    @test p4.xs[3] === p.xs[1]
    @test p4.xs[4] === p.xs[2]
    @test p4 == p+p

    # -
    m1 = e - t
    @test m1 isa Minus
    @test m1.p === e
    @test m1.q === t
    @test m1 == e-t
    m2 = -e
    @test m2 isa UniMinus
    @test m2.p === e
    @test m2 == -e

    # *
    t1 =  e * p
    @test t1 isa Times
    @test t1.xs[1] === e
    @test t1.xs[2] === p
    @test t1 == e*p
    t2 = t * e
    @test t2 isa Times
    @test t2.xs[1] === t.xs[1]
    @test t2.xs[2] === t.xs[2]
    @test t2.xs[3] === e
    @test t2 == t*e
    t3 = e * t
    @test t3 isa Times
    @test t3.xs[1] === e
    @test t3.xs[2] === t.xs[1]
    @test t3.xs[3] === t.xs[2]
    @test t3 == e*t
    t4 = t * t
    @test t4 isa Times
    @test t4.xs[1] === t.xs[1]
    @test t4.xs[2] === t.xs[2]
    @test t4.xs[3] === t.xs[1]
    @test t4.xs[4] === t.xs[2]
    @test t4 == t*t

    # adhoc *
    at1 = 3 * p
    @test at1 isa Times
    at2 = big(3) * p
    @test at2 isa Times
    @test at1 == at2
    at3 = p * 3
    @test at3 isa Times
    at4 = p * big(3)
    @test at4 isa Times
    @test at3 == at4

    # adhoc +
    ap1 = 3 + p
    @test ap1 isa Plus
    ap2 = big(3) + p
    @test ap2 isa Plus
    @test ap1 == ap2
    ap3 = p + 3
    @test ap3 isa Plus
    ap4 = p + big(3)
    @test ap4 isa Plus
    @test ap3 == ap4

    # adhoc -
    am1 = 3 - p
    @test am1 isa Minus
    am2 = big(3) - p
    @test am2 isa Minus
    @test am1 == am2
    am3 = p - 3
    @test am3 isa Minus
    am4 = p - big(3)
    @test am4 isa Minus
    @test am3 == am4

    # ^
    e1 = p^3
    @test e1 isa Exp
    @test e1.p === p
    @test e1.e == 3
    @test e1 == p^3

    h = Gen(:y)
    q = e1+t4*h
    @test gens(q) == [:x, :y]
    @test h == y
    @test q == e1+t4*h == ((1 + x)^3 + (x*(1 + x)*x*(1 + x)*y))
end

@testset "LazyPoly" begin
    F = LazyPolyRing(zz)
    r = Const(1) + Gen(:x)
    p = LazyPoly(F, r)
    @test parent(p) === F
    @test string(p) == "(1 + x)"
    for x in (gen(F, :x), F(:x))
        @test x isa LazyPoly{Int}
        @test x.p isa Gen
        @test x.p.g == :x
    end
    c1 = F(2)
    @test c1 isa LazyPoly{Int}
    @test c1.p isa Const{Int}
    @test c1.p.c === 2

    @test (p+c1).p isa Plus
    @test (p-c1).p isa Minus
    @test (-p).p isa UniMinus
    @test (p*c1).p isa Times
    @test (p^3).p isa Exp
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
    p = SLPoly(S, SLProgram{Int}())
    @test p isa SLPoly{Int,typeof(S)} <: MPolyElem{Int}
    @test parent(p) === S
    p = SLPoly(S)
    @test p isa SLPoly{Int,typeof(S)} <: MPolyElem{Int}
    @test parent(p) === S

    # copy
    q = SLPoly(S)
    # TODO: do smthg more interesting with q
    push!(constants(q), 3)
    push!(lines(q), Line(0))
    copy!(p, q)
    p2 = copy(q)
    for p1 in (p, p2)
        @test constants(p1) == constants(q) && constants(p1) !== constants(q)
        @test lines(p1) == lines(q) && lines(p1) !== lines(q)
    end
    S2 = SLPolyRing(zz, [:z, :t])
    @test_throws ArgumentError copy!(SLPoly(S2, SLProgram{Int}()), p)

    # building
    p = SLPoly(S)
    plines = lines(p)
    l1 = pushconst!(p.slprogram, 1)
    @test constants(p) == [1]
    @test l1 === SL.asconstant(1)

    # currently not supported anymore
    # l2 = pushconst!(p, 3)
    # @test constants(p) == [1, 3]
    # @test l2 === SL.asconstant(2)

    l3 = pushop!(p, SL.plus, l1, SL.input(1))
    @test l3 == Arg(UInt64(1))
    @test plines[1].x == 0x0340000018000001
    l4 = pushop!(p, SL.times, l3, SL.input(2))
    @test l4 == Arg(UInt64(2))
    @test plines[2].x ==0x0500000018000002
    pl = copy(plines)
    @test p === SL.pushfinalize!(p, l4)

    @test plines == [Line(0x0340000018000001), Line(0x0500000018000002)]
    SL.pushinit!(p)
    @test pl == plines
    SL.pushfinalize!(p, l4)
    @test plines == [Line(0x0340000018000001), Line(0x0500000018000002)]
    # p == (1+x)*y
    @test SL.evaluate!(Int[], p, [1, 2]) == 4
    @test SL.evaluate!(Int[], p, [0, 3]) == 3
    l5 = SL.pushinit!(p)
    l6 = SL.pushop!(p, SL.times, SL.input(1), SL.input(2)) # xy
    l7 = SL.pushop!(p, SL.exponentiate, l5, Arg(2)) # ((1+x)y)^2
    l8 = SL.pushop!(p, SL.minus, l6, l7) # xy - ((1+x)y)^2
    SL.pushfinalize!(p, l8)
    @test string(p) == "((xy) - ((1 + x)y)^2)"
    @test SL.evaluate!(Int[], p, [2, 3]) == -75
    @test SL.evaluate!(Int[], p, [-2, -1]) == 1

    @test SL.evaluate(p, [2, 3]) == -75
    @test SL.evaluate(p, [-2, -1]) == 1

    # compile!
    pf = SL.compile!(p)
    @test pf([2, 3]) == -75
    @test pf([-2, -1]) == 1
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
    @test q == -x1^2*y1^2-2*x1*y1^2+x1*y1-y1^2
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
    r = R()
    @test convert(R, convert(S, r)) == r

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
    p = S(x*y-16*y^2)
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

    for x = rand(Int64(0):Int(SL.cstmark-1), 100)
        if SL.isinput(Arg(x))
            @test SL.input(x) == Arg(x)
        else
            @test SL.input(x).x ⊻ SL.inputmark == x
        end
    end
end

@testset "SLProgram" begin
    p = SLProgram{Int}()
    @test p isa SLProgram{Int}
    @test isempty(p.cs)
    @test isempty(p.lines)
    @test !isassigned(p.f)

    # construction/evaluate/ninputs/aslazy
    p = SLProgram{Int}(1)
    @test evaluate(p, [10, 20]) == 10
    @test SL.ninputs(p) == 1
    @test SL.aslazy(p) == Gen(:x)
    p = SLProgram{Int}(3)
    @test evaluate(p, [10, "20", 'c']) == 'c'
    @test SL.ninputs(p) == 3
    @test SL.aslazy(p) == Gen(:z)

    p = SLProgram(Const(3))
    @test evaluate(p, [10, 20]) == 3
    @test SL.aslazy(p) == Const(3)
    p = SLProgram(Const('c'))
    @test evaluate(p, ["10", 20]) == 'c'
    @test SL.ninputs(p) == 0
    @test SL.aslazy(p) == Const('c')

    p = SLProgram{Int}(1)
    q = SLProgram(Const(6))
    r = SLProgram{Int}(2)
    x, y, z = Gen.([:x, :y, :z])

    # mutating ops
    @test p === SL.addeq!(p, q)
    @test evaluate(p, [3]) == 9
    @test SL.aslazy(p) == x+6

    @test p === SL.subeq!(p, r)
    @test evaluate(p, [3, 2]) == 7
    @test SL.aslazy(p) == x+6-y

    @test p === SL.subeq!(p)
    @test evaluate(p, [3, 2]) == -7
    @test SL.aslazy(p) == -(x+6-y)

    @test p === SL.muleq!(p, r)
    @test evaluate(p, [3, 2]) == -14
    @test SL.aslazy(p) == -(x+6-y)*y

    @test p === SL.expeq!(p, 3)
    @test evaluate(p, [3, 2]) == -2744
    @test SL.evaluates(p, [3, 2]) == [9, 7, -7, -14, -2744]
    @test SL.aslazy(p) == (-(x+6-y)*y)^3

    @test SL.ninputs(p) == 2

    p = SLProgram{UInt8}(1)
    q = SLProgram(Const(2))

    SL.addeq!(p, q)
    @test p.cs[1] === 0x2
    @test SL.aslazy(p) == x+2

    SL.muleq!(p, SLProgram(Const(3.0)))
    @test p.cs[2] === 0x3
    @test SL.aslazy(p) == (x+2)*3.0

    SL.subeq!(p, SLProgram(Const(big(4))))
    @test p.cs[3] === 0x4
    @test SL.aslazy(p) == (x+2)*3.0-big(4)

    @test_throws InexactError SL.addeq!(p, SLProgram(Const(1.2)))
    @assert length(p.cs) == 4 # p.cs was resized before append! failed
    pop!(p.cs) # set back consistent state
    @assert length(p.lines) == 4 # p.lines was resized before append! failed
    pop!(p.lines) # set back consistent state
    @test SL.aslazy(p) == (x+2)*3.0-big(4)

    p2 = SL.copy_oftype(p, Float64)
    @test p2 == p
    @test p2.cs == p.cs
    @test p2.lines == p.lines
    SL.addeq!(p2, SLProgram(Const(1.2)))
    @test p2.cs[4] == 1.2
    @test SL.aslazy(p2) == ((((x + 2.0)*3.0) - 4.0) + 1.2)

    p3 = copy(p)
    @test p3 == p
    @test p3.cs == p.cs
    @test p3.lines == p.lines
    @test_throws InexactError SL.addeq!(p3, SLProgram(Const(1.2)))

    # unary/binary ops
    p = SLProgram{BigInt}(1)
    q = SLProgram(Const(2))

    r = p+q
    @test SL.aslazy(r) == x+2
    @test SL.constantstype(r) === Signed

    r = r*SLProgram(Const(0x3))
    @test SL.aslazy(r) == (x+2)*3
    @test SL.constantstype(r) === Integer

    r = r-SLProgram(Const(1.2))
    @test SL.aslazy(r) == (x+2)*3-1.2
    @test SL.constantstype(r) === Real

    r = -r
    @test SL.aslazy(r) == -((x+2)*3-1.2)
    @test SL.constantstype(r) === Real

    r = r^3
    @test SL.aslazy(r) == (-((x+2)*3-1.2))^3
    @test SL.constantstype(r) === Real

    # conversion Lazy -> SLProgram
    @test SLProgram(x^2+y) isa SLProgram{Union{}}
    p = SL.muleq!(SLProgram(Const(2)), SLProgram{Int}(x^2+y))
    @test p isa SLProgram{Int}
    @test evaluate(p, [2, 3]) == 14
    @test SL.aslazy(p) == 2*(x^2 + y)
end
