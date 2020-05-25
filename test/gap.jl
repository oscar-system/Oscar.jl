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

@testset "GAPSLProgram compile!" begin
    x, y, z = slpgens(3)
    g = GAPSLProgram([ [2,3], [ 1, 2, 3, 1] ], 2 )
    slp = SL.compile!(g)
    @test g.slp[] === slp
    @test evaluate(slp, [x, y]) == x^2*y^3
    g = GAPSLProgram([ [2,3], [ 3, 1, 1, 4] ], 2 )
    slp = SL.compile!(g)
    @test evaluate(slp, [x, y]) == y^3*x^4
end
