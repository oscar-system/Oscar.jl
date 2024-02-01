@testset "Graphs" begin

    @testset "core functionality" begin
        g = Graph{Directed}(5)
        @test nvertices(g) == 5
        @test nedges(g) == 0
        add_edge!(g, 1, 2)
        @test nedges(g) == 1
        @test has_edge(g, 1, 2)
        rem_edge!(g, 1, 2)
        @test nedges(g) == 0
        @test !has_edge(g, 1, 2)
        @test add_vertex!(g)
        @test nvertices(g) == 6
        @test has_vertex(g, 6)
        rem_vertex!(g, 1)
        @test nvertices(g) == 5
        @test has_vertex(g, 1)
        @test !has_vertex(g, 6)
        @test add_vertices!(g, 5) == 5
        @test nvertices(g) == 10

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
        @test nvertices(egtriangle) == 3
        @test nedges(egtriangle) == 3
        @test nvertices(dgtriangle) == 3
        @test nedges(dgtriangle) == 3
        @test nvertices(egcube) == 8
        @test nedges(egcube) == 12
        @test nvertices(dgcube) == 6
        @test nedges(dgcube) == 12

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

        G = matrix(ZZ, 3, 3, [0,1,0,1,0,1,0,1,0])
        J = [2,3,1]
        H = G[J,J]
        b, I = Oscar._is_equal_up_to_permutation_with_permutation(G, H)
        @assert G[I,I] == H
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

    @testset "errors" begin
        g = Graph{Undirected}(1)
        @test !add_edge!(g,1,2)
    end
end
