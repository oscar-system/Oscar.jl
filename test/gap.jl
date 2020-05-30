@testset "GAPSLProgram construction" begin
    lines = GAPStraightLine[]

    pushline!(lines, [1, 2])
    @test lines[1] == [1, 2]
    @test_throws ArgumentError pushline!(lines, [1, 2, 3])

    pushline!(lines, ([2, 3], 1))
    @test lines[2] == ([2, 3], 1)
    @test_throws MethodError pushline!(lines, ([2, 3], 2, 3))

    pushline!(lines, [[4, 3], 2])
    @test lines[3] == ([4, 3], 2)
    @test_throws ArgumentError pushline!(lines, [[2, 3], 2, 3])
    @test_throws ArgumentError pushline!(lines, [[2, 3, 4], 2, 3])

    pushline!(lines, [[4, 3], [1, 2, 3, 4]])
    @test lines[4] == [[4, 3], [1, 2, 3, 4]]
    @test_throws ArgumentError pushline!(lines, [[2, 3], [2, 3, 4]])

    # return list not at the end
    @test_throws ArgumentError pushline!(lines, [1, 2])
    @test_throws ArgumentError pushline!(lines, ([2, 3], 1))
    @test_throws ArgumentError pushline!(lines, [[4, 3], 2])
    @test_throws ArgumentError pushline!(lines, [[4, 3], [1, 2, 3, 4]])

    @test_throws ArgumentError GAPSLProgram(lines)
    p = GAPSLProgram(lines, 9)
    @test p.lines == lines
    q = GAPSLProgram([[1, 2],
                      ([2, 3], 1),
                      [[4, 3], 2],
                      [[4, 3], [1, 2, 3, 4]]], 9)
    @test q.lines == lines

    r = GAPSLProgram(
        Any[ [[1, 2], 3],
             [[3, 2], 2],
             [1, 2, 2, 1] ] )
    @test r.ngens == 1

    @test_throws ArgumentError GAPSLProgram(
        Any[ [[1, 2], 3],
             [[3, 2], 2],
             [1, 2, 4, 1] ], 2)

    @test_throws ArgumentError GAPSLProgram(
        Any[ [[1, 2], 3],
             [[3, 2], 2],
             [[1, 2, 4, 1]] ], 2 )

    r = GAPSLProgram(
        Any[ [[1, 2], 3],
             [[1, 2, 4, 2], 2],
             [1, 2, 2, 1] ])
    @test r.ngens == 4

    r = GAPSLProgram(
        Any[ [[2, 2], 3],
             [[3, 2], 2],
             [1, 2, 2, 1] ])
    @test r.ngens == 2

    r = GAPSLProgram(
        Any[ [[2, 3], 3],
             [[3, 2], 2],
             [[1, 2, 2, 1], [4, 1]] ])
    @test r.ngens == 4
end

@testset "GAPSLProgram evaluate / compile!" begin
    for xy = (slpgens(3), Lazy[Gen(:x), Gen(:y)])
        x, y = xy
        res = empty(xy)

        g = GAPSLProgram([ [2,3], [ 1, 2, 3, 1] ], 2 )

        @test evaluate(g, xy) == x^2*y^3
        @test evaluate!(res, g, xy) == x^2*y^3
        slp = SL.compile!(g)
        @test g.slp[] === slp
        @test evaluate(slp, xy) == x^2*y^3

        g = GAPSLProgram([ [2,3], [ 3, 1, 1, 4] ], 2 )
        @test evaluate(g, xy) == y^3*x^4
        @test evaluate!(res, g, xy) == y^3*x^4
        slp = SL.compile!(g)
        @test evaluate(slp, xy) == y^3*x^4

        g = GAPSLProgram([ [2,3], [ [3, 1, 1, 4], [ 1, 2, 3, 1]] ], 2 )
        @test evaluate(g, xy) == [y^3*x^4, x^2*y^3]
        @test evaluate!(res, g, xy) == [y^3*x^4, x^2*y^3]
        slp = SL.compile!(g)
        @test evaluate(slp, xy) == [y^3*x^4, x^2*y^3]

        g = GAPSLProgram(Any[ [ [1,1,2,2], 2 ], [2,3,1,1] ] )
        @test evaluate(g, xy) == (x*y^2)^3*x
        @test evaluate!(res, g, xy) == (x*y^2)^3*x
        slp = SL.compile!(g)
        @test evaluate(slp, xy) == (x*y^2)^3*x

        g = GAPSLProgram(Any[ [ [1,-1,2,-2], 2 ], [2,-3,1,1] ] )
        t = (x^-Int(1)*y^-Int(2))^-Int(3)*x
        @test evaluate(g, xy) == t
        @test evaluate!(res, g, xy) == t
        slp = SL.compile!(g)
        @test evaluate(slp, xy) == t

        # assignments
        g = GAPSLProgram(Any[ [[2, 3], 3] ])
        @test evaluate(g, xy) == y^3
        slp = SL.compile!(g)
        @test evaluate(slp, xy) == y^3

        g = GAPSLProgram(Any[ [[2, 3], 3], [[1, 2, 3, 2], 2] ])
        @test evaluate(g, xy) == x^2*(y^3)^2
        slp = SL.compile!(g)
        @test evaluate(slp, xy) == x^2*(y^3)^2
    end
end
