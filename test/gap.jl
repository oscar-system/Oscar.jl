@testset "GAPSLProgram construction" begin
    p = GAPSLProgram()
    @test isempty(p.lines)

    @test p === pushline!(p, [1, 2])
    @test p.lines[1] == [1, 2]
    @test_throws ArgumentError pushline!(p, [1, 2, 3])

    @test p === pushline!(p, ([2, 3], 1))
    @test p.lines[2] == ([2, 3], 1)
    @test_throws MethodError pushline!(p, ([2, 3], 2, 3))

    @test p === pushline!(p, [[4, 3], 2])
    @test p.lines[3] == ([4, 3], 2)
    @test_throws ArgumentError pushline!(p, [[2, 3], 2, 3])
    @test_throws ArgumentError pushline!(p, [[2, 3, 4], 2, 3])


    @test p === pushline!(p, [[4, 3], [1, 2, 3, 4]])
    @test p.lines[4] == [[4, 3], [1, 2, 3, 4]]
    @test_throws ArgumentError pushline!(p, [[2, 3], [2, 3, 4]])

    # return list not at the end
    @test_throws ArgumentError pushline!(p, [1, 2])
    @test_throws ArgumentError pushline!(p, ([2, 3], 1))
    @test_throws ArgumentError pushline!(p, [[4, 3], 2])
    @test_throws ArgumentError pushline!(p, [[4, 3], [1, 2, 3, 4]])
end
