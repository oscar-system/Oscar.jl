using Test, StraightLinePrograms, AbstractAlgebra

using StraightLinePrograms: Const, ExpPoly, Gen, MinusPoly, PlusPoly, RecPoly,
    TimesPoly, UniMinusPoly

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
end
