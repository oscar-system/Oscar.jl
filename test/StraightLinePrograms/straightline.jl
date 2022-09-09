@testset "LazyRec" begin
    x, y, z = Gen.([:x, :y, :z])
    xyz = Any[x, y, z]

    # Const
    c = Const(1)
    @test c isa Const{Int}
    @test c isa LazyRec
    @test string(c) == "1"
    @test isempty(SLP.gens(c))
    @test c == Const(1) == Const(0x1)
    @test c != Const(2)
    @test SLP.evaluate(c, rand(3)) == 1
    @test SLP.evaluate(c, xyz) == 1

    # Gen
    g = Gen(:x)
    @test g isa Gen
    @test g isa LazyRec
    @test string(g) == "x"
    @test SLP.gens(g) == [:x]
    @test g == x
    @test g != y
    @test SLP.evaluate(g, [2, 3, 4]) == 2
    @test SLP.evaluate(g, xyz) == x
    @test SLP.evaluate(g, [Gen(:a), Gen(:b), Gen(:x)]) == Gen(:a)

    # Plus
    p = Plus(c, g)
    @test p isa Plus <: LazyRec
    @test p.xs[1] == c && p.xs[2] == g
    @test string(p) == "(1 + x)"
    @test SLP.gens(p) == [:x]
    @test p == 1+x == 0x1+x
    @test p != 2+x && p != 1+y
    @test SLP.evaluate(p, [2]) == 3
    @test SLP.evaluate(p, xyz) == 1 + x

    # Minus
    m = Minus(p, g)
    @test m isa Minus <: LazyRec
    @test string(m) == "((1 + x) - x)"
    @test SLP.gens(m) == [:x]
    @test m == (1+x)-x
    @test m != (1+x)+x && m != x-x && m != (1+x)-y
    @test SLP.evaluate(m, [2]) == 1
    @test SLP.evaluate(m, xyz) == (1+x)-x

    # UniMinus
    u = UniMinus(p)
    @test u isa UniMinus <: LazyRec
    @test string(u) == "(-(1 + x))"
    @test SLP.gens(u) == [:x]
    @test u == -(1+x)
    @test u != (1+x) && u != -(1+y)
    @test SLP.evaluate(u, [2]) == -3
    @test SLP.evaluate(u, xyz) == -(1 + x)

    # Times
    t = Times(g, p)
    @test t isa Times <: LazyRec
    @test string(t) == "(x*(1 + x))"
    @test SLP.gens(t) == [:x]
    @test t == x*(1+x)
    @test t != (1+x)*x && t != y*(1+x) && t != x*(1+y)
    @test SLP.evaluate(t, [2]) == 6
    @test SLP.evaluate(t, xyz) == x*(1+x)

    # Exp
    e = Exp(p, 3)
    @test e isa Exp <: LazyRec
    @test string(e) == "(1 + x)^3"
    @test SLP.gens(e) == [:x]
    @test e == (1+x)^3
    @test e != (1+x)^4 && e != (1+y)^3
    @test SLP.evaluate(e, [2]) == 27
    @test SLP.evaluate(e, xyz) == (1+x)^3

    # Call
    c = Call((x, y) -> 2x+3y, [x-y, y])
    @test c isa Call <: LazyRec
    @test SLP.gens(c) == [:x, :y]
    @test c == Call(c.f, [x-y, y])
    @test c != Call((x, y) -> 2x+3y, [x-y, y]) # not same function
    @test c != Call(c.f, [x-y, 2y])
    @test SLP.evaluate(c, [2, 3]) == 7

    # +
    p1 =  e + t
    @test p1 isa Plus
    @test p1.xs[1] === e
    @test p1.xs[2] === t
    @test p1 == e+t
    @test SLP.evaluate(p1, xyz) == (1+x)^3 + x*(1+x)

    p2 = p + e
    @test p2 isa Plus
    @test p2.xs[1] === p.xs[1]
    @test p2.xs[2] === p.xs[2]
    @test p2.xs[3] === e
    @test p2 == p+e
    @test SLP.evaluate(p2, xyz) == (1 + x + (1 + x)^3)

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
    @test SLP.gens(q) == [:x, :y]
    @test h == y
    @test q == e1+t4*h == ((1 + x)^3 + (x*(1 + x)*x*(1 + x)*y))
    @test SLP.evaluate(q, [2, 3]) == 135
    @test SLP.evaluate(q, xyz) == ((1 + x)^3 + (x*(1 + x)*x*(1 + x)*y))
end

@testset "Lazy" begin
    x, y, z = xyz = SLP.lazygens(3)

    @test xyz == SLP.gens(SLP.Lazy, 3)

    xs = Float64[2, 3, 4]

    @test SLP.evaluate(x, xs) == 2
    @test SLP.evaluate(y, xs) == 3
    @test SLP.evaluate(z, xs) == 4

    @test SLP.evaluate(x, xyz) == x
    @test SLP.evaluate(y, xyz) == y
    @test SLP.evaluate(z, xyz) == z

    # constant
    p = Lazy(1)
    @test SLP.evaluate(p, rand(Int, rand(0:20))) === 1
    p = Lazy([1, 2, 4])
    @test SLP.evaluate(p, rand(Int, rand(0:20))) == [1, 2, 4]

    p = 9 + 3*x*y^2 + ((y+z+3-x-3)*2)^-2 * 100
    @test SLP.evaluate(p, xs) == 64
    @test SLP.evaluate(p, xyz) == p

    a, b = SLP.lazygens([:a, :bc])
    @test string(a) == "a" && string(b) == "bc"

    q1 = SLP.compile(SLProgram, p)
    @test SLP.evaluate(q1, xs) == 64
    q2 = SLP.compile(p)
    @test q2 == q1
    @test SLP.evaluate(q2, xs) == 64
    q3 = SLProgram(p)
    @test q3 == q1
    @test SLP.evaluate(q3, xs) == 64

    x1, x2, x3 = SLP.gens(SLProgram, 3)
    @test SLProgram(x) == slpgen(1) == x1
    @test SLProgram(y) == slpgen(2) == x2
    @test SLProgram(z) == slpgen(3) == x3

    # call
    fun2 = (x, y) -> 2x+3y
    c = call(fun2, x-y, y)
    @test c isa Lazy
    @test SLP.gens(c) == [:x, :y]
    @test c == call(fun2, x-y, y)
    @test c != call(fun2, x-y, 2y)
    @test SLP.evaluate(c, [2, 3]) == 7

    c = call(fun2, 1, 3)
    @test isempty(SLP.gens(c))
    @test SLP.evaluate(c, []) == 11

    # SLP.evaluate: caching of results
    counter = 0
    c = call(_ -> (counter += 1; 0), x)
    p = c+c*c
    SLP.evaluate(p, [1])
    @test counter == 1
    # TODO: should add tests for every LazyRec subtypes
end

@testset "SL internals" begin
    @test SLP.showop == Dict(SLP.assign       => "->",
                            SLP.plus         => "+",
                            SLP.uniminus     => "-",
                            SLP.minus        => "-",
                            SLP.times        => "*",
                            SLP.divide       => "/",
                            SLP.exponentiate => "^",
                            SLP.keep         => "keep",
                            SLP.decision     => "&",
                            SLP.getindex_    => "[]")
    @test length(SLP.showop) == 10 # tests all keys are distinct
    for op in keys(SLP.showop)
        @test SLP.isassign(op) == (op == SLP.assign)
        @test SLP.istimes(op) == (op == SLP.times)
        # ...
        @test (op.x & 0x8000000000000000 != 0) ==
            SLP.isquasiunary(op) ==
            (op ∈ (SLP.uniminus, SLP.exponentiate, SLP.keep))
        @test SLP.isunary(op) == (op ∈ (SLP.uniminus, SLP.keep))
    end

    # pack & unpack
    ops = SLP.Op.(rand(UInt64(0):UInt64(0xff), 100) .<< 62)
    is = rand(UInt64(0):SLP.argmask, 100)
    js = rand(UInt64(0):SLP.argmask, 100)
    @test SLP.unpack.(SLP.pack.(ops, Arg.(is), Arg.(js))) == tuple.(ops, Arg.(is), Arg.(js))

    for x = rand(Int64(0):Int(SLP.cstmark-1), 100)
        if SLP.isinput(Arg(x))
            @test SLP.input(x) == Arg(x)
        else
            @test SLP.input(x).x ⊻ SLP.inputmark == x
        end
    end
end

@testset "Arg" begin
    @test_throws InexactError Arg(-1)
    @test_throws ArgumentError Arg(typemax(Int))
    @test_throws ArgumentError Arg(1 + SLP.argmask % Int)
    a = Arg(SLP.argmask)
    @test a.x == SLP.argmask

    @test_throws ArgumentError SLP.intarg(SLP.payloadmask % Int)
    x = (SLP.payloadmask ⊻ SLP.negbit) % Int
    @test SLP.getint(SLP.intarg(x)) == x
    @test_throws ArgumentError SLP.intarg(x+1)
    @test SLP.getint(SLP.intarg(-x)) == -x
    @test SLP.getint(SLP.intarg(-x-1)) == -x-1
    @test_throws ArgumentError SLP.intarg(-x-2)
end

@testset "SLProgram" begin
    x, y, z = Gen.([:x, :y, :z])
    xyz = Any[x, y, z]

    p = SLProgram()
    @test p isa SLProgram{Union{}}
    @test isempty(p.cs)
    @test isempty(p.lines)
    @test !isassigned(p.f)

    p = SLProgram{Int}()
    @test p isa SLProgram{Int}
    @test isempty(p.cs)
    @test isempty(p.lines)
    @test !isassigned(p.f)

    # construction/SLP.evaluate/ninputs/aslazy
    p = SLProgram{Int}(1)
    @test SLP.evaluate(p, [10, 20]) == 10
    @test SLP.ninputs(p) == 1
    @test SLP.aslazyrec(p) == Gen(:x)
    p = SLProgram(3)
    @test SLP.evaluate(p, [10, "20", 'c']) == 'c'
    @test SLP.ninputs(p) == 3
    @test SLP.aslazyrec(p) == Gen(:z)

    p = SLProgram(Const(3))
    @test SLP.evaluate(p, [10, 20]) == 3
    @test SLP.aslazyrec(p) == Const(3)
    p = SLProgram(Const('c'))
    @test SLP.evaluate(p, ["10", 20]) == 'c'
    @test SLP.ninputs(p) == 0
    @test SLP.aslazyrec(p) == Const('c')

    # exponent
    p = SLProgram(x*y^-2)
    @test SLP.aslazyrec(p) == x*y^Int(-2)
    e = (SLP.payloadmask ⊻ SLP.negbit) % Int
    @test x^e == SLP.aslazyrec(SLProgram(x^e))
    @test_throws ArgumentError SLProgram(x^(e+1))
    e = -e-1
    @test x^e == SLP.aslazyrec(SLProgram(x^e))
    @test_throws ArgumentError SLProgram(x^(e-1))
    p = SLProgram(x^2*y^(Int(-3)))
    @test SLP.evaluate(p, [sqrt(2), 4^(-1/3)]) === 8.0
    @test SLP.evaluate(p^-1, [sqrt(2), 4^(-1/3)]) == 0.125
    p = SLProgram(x^0)
    @test SLP.evaluate(p, [2]) == 1
    @test SLP.evaluate(p, xyz) == x^0

    # assign
    p = SLProgram{Int}()
    k = SLP.pushop!(p, SLP.plus, SLP.input(1), SLP.input(2))
    k = SLP.pushop!(p, SLP.times, k, SLP.input(2))
    @assert length(p.lines) == 2
    k = SLP.pushop!(p, SLP.assign, k, Arg(1))
    @test k == Arg(1)
    SLP.pushfinalize!(p, k)
    @test SLP.evaluate(p, LazyRec[x, y]) == (x+y)*y
    k = SLP.pushop!(p, SLP.exponentiate, k, Arg(2))
    @test SLP.evaluate(p, LazyRec[x, y]) == (x+y)*y
    SLP.pushfinalize!(p, k)
    @test SLP.evaluate(p, LazyRec[x, y]) == ((x+y)*y)^2
    SLP.pushfinalize!(p, Arg(2))
    @test SLP.evaluate(p, LazyRec[x, y]) == ((x+y)*y)
    k = SLP.pushop!(p, SLP.assign, k, Arg(2))
    @test k == Arg(2)
    @test SLP.evaluate(p, LazyRec[x, y]) == ((x+y)*y)^2

    # nsteps
    @test nsteps(p) == 5

    # permute_inputs!
    x1, x2, x3 = slpgens(3)
    p = x1*x2^2+x3^3
    for perm = ([3, 1, 2], Perm([3, 1, 2]))
        q = copy(p)
        SLP.permute_inputs!(q, perm)
        @test q == x3*x1^2+x2^3
    end

    # nsteps again
    @test nsteps(p) == 4
    @test nsteps(x1) == nsteps(x3) == 0

    p = SLProgram()
    i = SLP.pushint!(p, 2)
    j = SLP.pushint!(p, 4)
    k = SLP.pushop!(p, SLP.plus, i, j)
    SLP.pushfinalize!(p, k)
    @test nsteps(p) == 1

    # mutating ops
    p = SLProgram{Int}(1)
    q = SLProgram(Const(6))
    r = SLProgram(2)

    @test p === SLP.addeq!(p, q)
    @test SLP.evaluate(p, [3]) == 9
    @test SLP.aslazyrec(p) == x+6

    @test p === SLP.subeq!(p, r)
    @test SLP.evaluate(p, [3, 2]) == 7
    @test SLP.aslazyrec(p) == x+6-y

    @test p === SLP.subeq!(p)
    @test SLP.evaluate(p, [3, 2]) == -7
    @test SLP.aslazyrec(p) == -(x+6-y)

    @test p === SLP.muleq!(p, r)
    @test SLP.evaluate(p, [3, 2]) == -14
    @test SLP.aslazyrec(p) == -(x+6-y)*y

    @test p === SLP.expeq!(p, 3)
    @test SLP.evaluate(p, [3, 2]) == -2744
    @test SLP.evaluates(p, [3, 2]) == [9, 7, -7, -14, -2744]
    @test SLP.aslazyrec(p) == (-(x+6-y)*y)^3

    @test SLP.ninputs(p) == 2

    p = SLProgram{UInt8}(1)
    q = SLProgram(Const(2))

    SLP.addeq!(p, q)
    @test p.cs[1] === 0x2
    @test SLP.aslazyrec(p) == x+2

    SLP.muleq!(p, SLProgram(Const(3.0)))
    @test p.cs[2] === 0x3
    @test SLP.aslazyrec(p) == (x+2)*3.0

    SLP.subeq!(p, SLProgram(Const(big(4))))
    @test p.cs[3] === 0x4
    @test SLP.aslazyrec(p) == (x+2)*3.0-big(4)

    @test_throws InexactError SLP.addeq!(p, SLProgram(Const(1.2)))
    @assert length(p.cs) == 4 # p.cs was resized before append! failed
    pop!(p.cs) # set back consistent state
    @assert length(p.lines) == 3 # p.lines was *not* resized before append! failed
    @test SLP.aslazyrec(p) == (x+2)*3.0-big(4)

    p2 = SLP.copy_oftype(p, Float64)
    @test p2 == p
    @test p2.cs == p.cs
    @test p2.lines == p.lines
    SLP.addeq!(p2, SLProgram(Const(1.2)))
    @test p2.cs[4] == 1.2
    @test SLP.aslazyrec(p2) == ((((x + 2.0)*3.0) - 4.0) + 1.2)

    p3 = copy(p)
    @test p3 == p
    @test p3.cs == p.cs
    @test p3.lines == p.lines
    @test_throws InexactError SLP.addeq!(p3, SLProgram(Const(1.2)))

    # unary/binary ops
    p = SLProgram{BigInt}(1)
    p2 = SLProgram(1)
    q = SLProgram(Const(2))

    r = p+q
    @test SLP.aslazyrec(r) == x+2
    @test SLP.constantstype(r) === Signed

    r2 = p2+q
    @test SLP.aslazyrec(r) == x+2
    @test SLP.constantstype(r2) === Int

    r = r*SLProgram(Const(0x3))
    @test SLP.aslazyrec(r) == (x+2)*3
    @test SLP.constantstype(r) === Integer

    r2 = r2*SLProgram(Const(0x3))
    @test SLP.aslazyrec(r2) == (x+2)*3
    @test SLP.constantstype(r2) === Integer

    r = r-SLProgram(Const(1.2))
    @test SLP.aslazyrec(r) == (x+2)*3-1.2
    @test SLP.constantstype(r) === Real

    r = -r
    @test SLP.aslazyrec(r) == -((x+2)*3-1.2)
    @test SLP.constantstype(r) === Real

    r = r^3
    @test SLP.aslazyrec(r) == (-((x+2)*3-1.2))^3
    @test SLP.constantstype(r) === Real

    @testset "adhoc" begin
        r = p+q
        @test typeof(r) == SLProgram{Signed}
        @assert SLP.evaluate(r, xyz) == x+2

        @test SLP.evaluate(2+r, xyz) == 2+(x+2)
        @test SLP.evaluate(r+big(4), xyz) == (x+2) + big(4)

        @test SLP.evaluate(2*r, xyz) == 2*(x+2)
        @test SLP.evaluate(r*1.3, xyz) == (x+2) * 1.3

        @test SLP.evaluate(2 - r, xyz) == 2 - (x+2)
        @test SLP.evaluate(r - 0x12, xyz) == (x+2) - 0x12
    end

    # conversion LazyRec -> SLProgram
    @test SLProgram(x^2+y) isa SLProgram{Union{}}
    p = SLP.muleq!(SLProgram(Const(2)), SLProgram{Int}(x^2+y))
    @test p isa SLProgram{Int}
    @test SLP.evaluate(p, [2, 3]) == 14
    @test SLP.aslazyrec(p) == 2*(x^2 + y)

    l = SLP.test(x^2 * y, 2) & SLP.test(y^2 * x, 3)
    p = SLProgram(l)
    P = AbstractAlgebra.SymmetricGroup(4)
    p1, p2 = P("(1,4,3)"), P("(1, 3)")
    @test SLP.evaluate(p, [p1, p2])
    @test !SLP.evaluate(p, [p2, p1])
    @test SLP.aslazyrec(p) == l

    # multiple return
    p = SLProgram{Int}()
    inputs = Any[x, y, z]

    @test SLP.evaluate(p, inputs) == []
    k1 = SLP.pushop!(p, SLP.plus, SLP.input(1), SLP.input(2))
    @test SLP.evaluate(p, inputs) == [x+y]
    k2 = SLP.pushconst!(p, 3)
    k3 = SLP.pushop!(p, SLP.times, k1, k2)
    @test SLP.evaluate(p, inputs) == [x+y, (x+y)*3]
    SLP.pushop!(p, SLP.assign, k3, k1)
    @test SLP.evaluate(p, inputs) == [(x+y)*3, (x+y)*3]
    SLP.pushfinalize!(p, k3)
    @test SLP.evaluate(p, inputs) == (x+y)*3
    SLP.setmultireturn!(p)
    @test SLP.evaluate(p, inputs) == [(x+y)*3, (x+y)*3]

    # multiple return & list
    X, Y = slpgens(2)
    pl = SLP.list([X*Y, X+1-Y])
    @test SLP.evaluate(pl, inputs) == [(x*y), (x+1-y)]
    pl = SLP.list([X, Y, X+Y]) # first elements don't add a "step"/"line"
    @test SLP.evaluate(pl, inputs) == [x, y, x+y]

    @test SLP.evaluate(X^2*Y+Y^2, [X, Y]) == X^2*Y+Y^2

    # keep
    @test_throws ArgumentError SLP.pushop!(p, SLP.keep, Arg(3))
    # test we are still in valid state:
    @test SLP.evaluate(p, inputs) == [(x+y)*3, (x+y)*3]
    SLP.pushop!(p, SLP.keep, Arg(2))
    @test SLP.evaluate(p, inputs) == [(x+y)*3, (x+y)*3]
    SLP.pushop!(p, SLP.keep, Arg(1))
    @test SLP.evaluate(p, inputs) == [(x+y)*3]
    SLP.pushop!(p, SLP.keep, Arg(0))
    @test SLP.evaluate(p, inputs) == []

    # integers
    p = SLProgram()
    i = SLP.pushint!(p, 123)
    j = SLP.pushint!(p, -4)

    k = SLP.pushop!(p, SLP.plus, i, j)
    SLP.pushfinalize!(p, k)
    @test SLP.evaluate(p, inputs) == 119

    k = SLP.pushop!(p, SLP.times, k, SLP.input(1))
    SLP.pushfinalize!(p, k)
    @test SLP.evaluate(p, inputs) == 119 * x

    l = SLP.pushint!(p, -2)
    SLP.pushfinalize!(p, l)
    @test SLP.evaluate(p, inputs) == -2

    m = SLP.pushop!(p, SLP.minus, k, l)
    SLP.pushfinalize!(p, m)
    @test SLP.evaluate(p, inputs) == 119 * x - (-2)

    @testset "bug with ints" begin
        u = SLProgram()
        i = SLP.pushint!(u, 1)
        j = SLP.pushint!(u, 2)
        k = SLP.pushop!(u, SLP.plus, i, j)
        SLP.pushfinalize!(u, k)
        v = SLProgram()
        i = SLP.pushint!(v, 3)
        j = SLP.pushint!(v, 4)
        k = SLP.pushop!(v, SLP.plus, i, j)
        SLP.pushfinalize!(v, k)
        w = u + v
        @test SLP.evaluate(w, []) == 10
    end

    # 3-args SLP.evaluate
    x, y = slpgens(2)
    p = 3*x*y^3
    s = SLP.evaluate(p, ["a", "b"], string)
    @test s == "3abbb"
    s = SLP.evaluate!(String[], p, ["a", "b"], string)
    @test s == "3abbb"

    p = [1, 3, 2]*x + (1, 2)
    @test SLP.evaluate(p, [2], length) == 8
end

@testset "SL Decision" begin
    x, y = SLP.lazyrecgens(2)
    a, b = ab = SLP.gens(SLP.Lazy, 2)

    p = SLProgram()
    pushop!(p, SLP.decision, SLP.input(1), SLP.pushint!(p, 3))
    c = pushop!(p, SLP.times, SLP.input(1), SLP.input(2))
    pushop!(p, SLP.decision, c, SLP.pushint!(p, 2))
    SLP.setdecision!(p)

    l = SLP.test(x, 3) & SLP.test(x*y, 2)
    f = SLP.test(a, 3) & SLP.test(a*b, 2)
    @test SLP.evaluate(p, Any[x, y]) == l
    @test SLP.evaluate(l, Any[x, y]) == l
    @test SLP.evaluate(l, SLP.gens(SLProgram, 2)) == p
    @test SLP.evaluate(p, ab) == f
    @test SLP.evaluate(f, slpgens(2)) == p

    S = SymmetricGroup(4)
    for (x, y) in eachcol(rand(S, 2, 200))
        res = order(x) == 3 && order(x*y) == 2
        @test SLP.evaluate(p, [x, y]) == res
        @test SLP.evaluate(l, [x, y]) == res
    end
end

@testset "SL lists" begin
    x, y = SLP.lazygens(2)
    X, Y = slpgens(2)

    q = SLP.list([x*y^2, x+1-y])
    @test SLP.evaluate(q, [x, y]) == q
    @test SLP.evaluate(q, Any[x, y]) == [x*y^2, x+1-y]
    @test SLP.evaluate(q, [2, 3]) == [18, 0]
    @test SLP.evaluate(SLP.evaluate(q, SLProgram[X, Y]), [x, y]) == q

    @test SLP.evaluate(SLP.list([q, q]), [x, y]) == SLP.list([q, q])
    # TODO: list of list of SLProgram hits an assertion error, so
    # handle this case more gracefully

    # test 3+ elements
    r = SLP.list([x, y+1, x+y, y-x])
    @test SLP.evaluate(r, [2, 3]) == [2, 4, 5, 1]
    @test SLP.evaluate(r, [x, y]) == r
end

@testset "SL compose" begin
    x, y = SLP.lazygens(2)
    X, Y = slpgens(2)

    q = SLP.compose(x - y, SLP.list([y, x]))
    p = SLP.compose(x - y, SLP.list([y, x]), flatten=false)
    r = SLP.compose(X - Y, SLP.list([Y, X]))

    for s = (q, p, r)
        @test SLP.evaluate(s, [2, 3]) == 1
        @test SLP.evaluate(s, [3.0, 1.0]) == -2
    end
    @test SLP.evaluate(r, [x, y]) == q

    q = SLP.compose(SLP.list([x+y, x-y]), SLP.list([y-x, y+x]))
    r = SLP.compose(SLP.list([X+Y, X-Y]), SLP.list([Y-X, Y+X]))

    for s = (q, r)
        @test SLP.evaluate(s, [2, 3]) == [6, -4]
        @test SLP.evaluate(s, [3.0, 1.0]) == [2, -6]
    end
    @test SLP.evaluate(r, [x, y]) == q

    p = SLP.compose(2.0*X+1.0, SLP.list([3*X]))
    @test p isa SLProgram{Real}
    @test SLP.evaluate(p, [x]) == 2*(3*x)+1
end

@testset "SL getindex" begin
    x, y = SLP.lazygens(2)
    X, Y = slpgens(2)

    p = y[x+1]
    @test p.x isa SLP.Getindex
    c = compile(SLProgram, p)
    @test c isa SLProgram{Int}
    q = Y[X+1]
    @test q == c
    @test q isa SLProgram{Int}

    for r = (p, q)
        @test SLP.evaluate(p, Any[2, [4, 5, 6]]) == 6
    end

    # adhoc
    p = y[2] + x[big(1)]
    @test SLP.evaluate(p, [[1, 2, 3], [10, 11, 12]]) == 12
    p = y[x[3]]
    @test SLP.evaluate(p, Any[[1, 2, 4], 1:4]) == 4

    p = SLP.slpcst([10, 20, 30])[2]
    @test SLP.evaluate(p, []) == 20
    p = SLP.list([X, Y, SLP.slpcst(30)])
    @test p isa SLProgram
    @test p[2] isa SLProgram
    @test p[3] isa SLProgram
    @test SLP.evaluate(p[2], [10, 20]) == 20
    @test SLP.evaluate(p[0x3], [10, 20]) == 30

    # multi-indices
    p = y[x, 2, 1]
    @test SLP.evaluate(p, [3, reshape(1:27, 3, 3, 3)]) == 6

    # integer & integer-array indexing
    p = list([x, y, y-x, y+x])

    @test p[1] == x
    @test p[0x2] == y
    @test p[big(4)] == y+x

    @test p[[1, 3]] == list([x, y-x])
    @test p[Any[1, 3]] == list([x, y-x])
    @test p[[4]] == list([y+x])
    @test p[Number[4]] == list([y+x])

    p = Lazy(Vector{Char})[['a']]
    e = SLP.evaluate(p, [])
    @test e isa Vector{Vector{Char}}
    @test e == [['a']]
end
