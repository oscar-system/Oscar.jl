using Test, StraightLinePrograms, AbstractAlgebra

using StraightLinePrograms: Const, Gen, PlusPoly, RecPoly

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
end
