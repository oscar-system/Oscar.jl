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
end
