@testset "Graphs" begin

    @testset "core functionality" begin
        g = Graph{Directed}(5)
        @test nv(g) == 5
        @test ne(g) == 0
        add_edge!(g, 1, 2)
        @test ne(g) == 1
        @test has_edge(g, 1, 2)
        rem_edge!(g, 1, 2)
        @test ne(g) == 0
        @test !has_edge(g, 1, 2)
        v = add_vertex!(g)
        @test v == 6
        @test nv(g) == 6
        @test has_vertex(g, 6)
        rem_vertex!(g, 1)
        @test nv(g) == 5
        @test has_vertex(g, 1)
        @test !has_vertex(g, 6)
    end

    triangle = simplex(2)
    c = cube(3)
    egtriangle = edgegraph(triangle)
    dgtriangle = dualgraph(triangle)
    egcube = edgegraph(c)
    dgcube = dualgraph(c)
    
    @testset "graphs from polytopes" begin
        @test nv(egtriangle) == 3
        @test ne(egtriangle) == 3
        @test nv(dgtriangle) == 3
        @test ne(dgtriangle) == 3
        @test nv(egcube) == 8
        @test ne(egcube) == 12
        @test nv(dgcube) == 6
        @test ne(dgcube) == 12
    end

end
