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

        g = Graph{Directed}(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        @test signed_incidence_matrix(g) == Matrix([-1 0 1; 1 -1 0; 0 1 -1; 0 0 0])
    end

    triangle = simplex(2)
    c = cube(3)
    cr = cross_polytope(3)
    egtriangle = edgegraph(triangle)
    dgtriangle = dualgraph(triangle)
    egcube = edgegraph(c)
    dgcube = dualgraph(c)
    egcr = edgegraph(cr)
    
    @testset "graphs from polytopes" begin
        @test nv(egtriangle) == 3
        @test ne(egtriangle) == 3
        @test nv(dgtriangle) == 3
        @test ne(dgtriangle) == 3
        @test nv(egcube) == 8
        @test ne(egcube) == 12
        @test nv(dgcube) == 6
        @test ne(dgcube) == 12

        @test is_isomorphic(dgtriangle, egtriangle)

        @test is_isomorphic(egcr, dgcube)
        @test !is_isomorphic(egcr, egcube)

        @test incidence_matrix(egtriangle) == IncidenceMatrix([[1,2],[1,3],[2,3]])
    end

    @testset "isomorphic" begin
        g = Graph{Directed}(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)

        gg = Graph{Directed}(5)
        add_edge!(gg, 3, 5)
        add_edge!(gg, 5, 4)
        @test !is_isomorphic(g, gg)

        add_edge!(gg, 3, 4)
        @test !is_isomorphic(g, gg)

        rem_edge!(gg, 3, 4)
        add_edge!(gg, 4, 3)
        @test is_isomorphic(g, gg)
    end

    @testset "connectivity" begin
        g = Graph{Directed}(5)
        @test !is_weakly_connected(g)
        @test !is_strongly_connected(g)
        @test length(weakly_connected_components(g)) == 5
        @test length(strongly_connected_components(g)) == 5

        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        add_edge!(g, 3, 4)
        add_edge!(g, 5, 3)
        @test is_weakly_connected(g)
        @test !is_strongly_connected(g)
        @test length(weakly_connected_components(g)) == 1
        @test length(strongly_connected_components(g)) == 3

        add_edge!(g, 4, 5)
        @test is_weakly_connected(g)
        @test is_strongly_connected(g)
        @test length(weakly_connected_components(g)) == 1
        @test length(strongly_connected_components(g)) == 1
        @test diameter(g) == 4

        g = Graph{Undirected}(5)
        @test !is_connected(g)
        @test length(connected_components(g)) == 5

        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        add_edge!(g, 4, 5)
        @test !is_connected(g)
        @test length(connected_components(g)) == 2

        add_edge!(g, 3, 5)
        @test is_connected(g)
        @test length(connected_components(g)) == 1
        @test diameter(g) == 3
    end
end
