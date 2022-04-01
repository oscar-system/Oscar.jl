@testset "GAPSLProgram construction" begin
    slines = GAPStraightLine[]

    pushline!(slines, [1, 2])
    @test slines[1] == [1, 2]
    @test_throws ArgumentError pushline!(slines, [1, 2, 3])

    pushline!(slines, ([2, 3], 1))
    @test slines[2] == ([2, 3], 1)
    @test_throws MethodError pushline!(slines, ([2, 3], 2, 3))

    pushline!(slines, [[4, 3], 2])
    @test slines[3] == ([4, 3], 2)
    @test_throws ArgumentError pushline!(slines, [[2, 3], 2, 3])
    @test_throws ArgumentError pushline!(slines, [[2, 3, 4], 2, 3])

    pushline!(slines, [[4, 3], [1, 2, 3, 4]])
    @test slines[4] == [[4, 3], [1, 2, 3, 4]]
    @test_throws ArgumentError pushline!(slines, [[2, 3], [2, 3, 4]])

    # return list not at the end
    @test_throws ArgumentError pushline!(slines, [1, 2])
    @test_throws ArgumentError pushline!(slines, ([2, 3], 1))
    @test_throws ArgumentError pushline!(slines, [[4, 3], 2])
    @test_throws ArgumentError pushline!(slines, [[4, 3], [1, 2, 3, 4]])

    @test_throws ArgumentError GAPSLProgram(slines)
    p = GAPSLProgram(slines, 9)
    @test p.lines == slines
    q = GAPSLProgram([[1, 2],
                      ([2, 3], 1),
                      [[4, 3], 2],
                      [[4, 3], [1, 2, 3, 4]]], 9)
    @test q.lines == slines

    r = GAPSLProgram(
        [ [[1, 2], 3],
          [[3, 2], 2],
          [1, 2, 2, 1] ] )
    @test r.ngens == 1

    @test_throws ArgumentError GAPSLProgram(
        [ [[1, 2], 3],
          [[3, 2], 2],
          [1, 2, 4, 1] ], 2)

    @test_throws ArgumentError GAPSLProgram(
        [ [[1, 2], 3],
          [[3, 2], 2],
          [[1, 2, 4, 1]] ], 2 )

    r = GAPSLProgram(
        [ [[1, 2], 3],
          [[1, 2, 4, 2], 2],
          [1, 2, 2, 1] ])
    @test r.ngens == 4

    r = GAPSLProgram(
        [ [[2, 2], 3],
          [[3, 2], 2],
          [1, 2, 2, 1] ])
    @test r.ngens == 2

    r = GAPSLProgram(
        [ [[2, 3], 3],
          [[3, 2], 2],
          [[1, 2, 2, 1], [4, 1]] ])
    @test r.ngens == 4
end

@testset "GAPSLDecision construction" begin
    d = GAPSLDecision([
        [ [ 1, 1, 2, 1 ], 3 ],
        (1, 2), # "Order" line
        (2, 3),
        (3, 5) ])
    @test d.ngens == 2

    d = GAPSLDecision([
        [ [ 1, 1, 2, 1 ], 3 ],
        ["Order", 1, 2],
        ["Order", 2, 3],
        ["Order", 3, 5] ])
    @test d.ngens == 2

    d = GAPSLDecision([
        [ 1, 1, 2, 1 ],
        (1, 2), # "Order" line
        (2, 3),
        (3, 5) ], 2)
    @test d.ngens == 2

    # ngens must be specified
    @test_throws ArgumentError GAPSLDecision([ [1, 1], (1, 1) ])
    @test_throws ArgumentError GAPSLDecision([ [1, 1], ["Order", 1, 1] ])
    # can't have return line
    @test_throws ArgumentError GAPSLDecision([ [[1, 1], 2], [[1, 1]] ])
end

@testset "GAPSLProgram evaluate / compile!" begin
    for xy = (slpgens(3), LazyRec[Gen(:x), Gen(:y)])
        x, y = xy
        res = empty(xy)

        g = GAPSLProgram([ [2,3], [ 1, 2, 3, 1] ], 2 )

        @test SLP.evaluate(g, xy) == x^2*y^3
        @test SLP.evaluate!(res, g, xy) == x^2*y^3
        slp0 = SLP.compile(g)
        @test slp0 isa SLProgram
        @test SLP.evaluate(slp0, xy) == x^2*y^3
        slp = SLP.compile(SLProgram, g)
        @test slp isa SLProgram
        @test slp == slp0
        @test SLP.evaluate(slp, xy) == x^2*y^3
        slp = SLP.compile!(g)
        @test slp == slp0
        @test g.slp[] === slp
        @test SLP.evaluate(slp, xy) == x^2*y^3

        g = GAPSLProgram([ [2,3], [ 3, 1, 1, 4] ], 2 )
        @test SLP.evaluate(g, xy) == y^3*x^4
        @test SLP.evaluate!(res, g, xy) == y^3*x^4
        slp = SLP.compile!(g)
        @test SLP.evaluate(slp, xy) == y^3*x^4

        g = GAPSLProgram([ [2,3], [ [3, 1, 1, 4], [ 1, 2, 3, 1]] ], 2 )
        @test SLP.evaluate(g, xy) == [y^3*x^4, x^2*y^3]
        @test SLP.evaluate!(res, g, xy) == [y^3*x^4, x^2*y^3]
        slp = SLP.compile!(g)
        # convert to Vector{Any} as otherwise we don't get a vector
        # but an SLP returning a list
        @test SLP.evaluate(slp, convert(Vector{Any}, xy)) == [y^3*x^4, x^2*y^3]

        g = GAPSLProgram([ [ [1,1,2,2], 2 ], [2,3,1,1] ] )
        @test SLP.evaluate(g, xy) == (x*y^2)^3*x
        @test SLP.evaluate!(res, g, xy) == (x*y^2)^3*x
        slp = SLP.compile!(g)
        @test SLP.evaluate(slp, xy) == (x*y^2)^3*x

        g = GAPSLProgram([ [ [1,-1,2,-2], 2 ], [2,-3,1,1] ] )
        t = (x^-Int(1)*y^-Int(2))^-Int(3)*x
        @test SLP.evaluate(g, xy) == t
        @test SLP.evaluate!(res, g, xy) == t
        slp = SLP.compile!(g)
        @test SLP.evaluate(slp, xy) == t

        # assignments
        g = GAPSLProgram([ [[2, 3], 3] ])
        @test SLP.evaluate(g, xy) == y^3
        slp = SLP.compile!(g)
        @test SLP.evaluate(slp, xy) == y^3

        g = GAPSLProgram([ [[2, 3], 3], [[1, 2, 3, 2], 2] ])
        @test SLP.evaluate(g, xy) == x^2*(y^3)^2
        slp = SLP.compile!(g)
        @test SLP.evaluate(slp, xy) == x^2*(y^3)^2
    end
end

@testset "GAPSLDecision evaluate / compile!" begin
    p = perm"(1, 2)(3, 4)"
    q = Perm([1, 4, 2, 3])
    pq = [p, q]

    for (i, j, k, r) in [(2, 3, 5, false),
                         (2, 3, 3, true),
                         (2, 2, 3, false),
                         (1, 3, 3, false)]
        gd = GAPSLDecision([
            [ [ 1, 1, 2, 1 ], 3 ],
            ["Order", 1, i],
            ["Order", 2, j],
            ["Order", 3, k] ])

        @test SLP.evaluate(gd, pq) == r

        gdc = SLP.compile!(gd)
        @test SLP.evaluate(gdc, pq) == r

        sld = SLProgram()
        m = SLP.pushop!(sld, SLP.times, SLP.input(1), SLP.input(2))
        SLP.pushop!(sld, SLP.decision, SLP.input(1), SLP.pushint!(sld, i))
        SLP.pushop!(sld, SLP.decision, SLP.input(2), SLP.pushint!(sld, j))
        dec = SLP.pushop!(sld, SLP.decision, m, SLP.pushint!(sld, k))
        SLP.pushfinalize!(sld, dec)
        @test SLP.evaluate(sld, pq) == r
    end
end
