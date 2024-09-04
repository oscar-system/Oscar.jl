@testset "Graphs" begin

    @testset "core functionality" begin
        g = Graph{Directed}(5)
        @test n_vertices(g) == 5
        @test n_edges(g) == 0
        add_edge!(g, 1, 2)
        @test n_edges(g) == 1
        @test has_edge(g, 1, 2)
        rem_edge!(g, 1, 2)
        @test n_edges(g) == 0
        @test !has_edge(g, 1, 2)
        @test add_vertex!(g)
        @test n_vertices(g) == 6
        @test has_vertex(g, 6)
        rem_vertex!(g, 1)
        @test n_vertices(g) == 5
        @test has_vertex(g, 1)
        @test !has_vertex(g, 6)
        @test add_vertices!(g, 5) == 5
        @test n_vertices(g) == 10
        @test rem_vertices!(g, [2, 4, 6, 11])
        @test n_vertices(g) == 7

        g = Graph{Directed}(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        @test signed_incidence_matrix(g) == Matrix([-1 0 1; 1 -1 0; 0 1 -1; 0 0 0])
    end

    triangle = simplex(2)
    c = cube(3)
    cr = cross_polytope(3)
    pos = convex_hull([0 0 0; 1 0 0], [0 1 0; 0 0 1])
    pl = convex_hull([0 0 0; 1 0 0], nothing, [0 1 1])
    egtriangle = vertex_edge_graph(triangle)
    dgtriangle = dual_graph(triangle)
    egcube = vertex_edge_graph(c)
    dgcube = dual_graph(c)
    egcr = vertex_edge_graph(cr)
    egpos = vertex_edge_graph(pos)
    egpl = vertex_edge_graph(pl)
    egplc = vertex_edge_graph(pl, modulo_lineality=true)
    
    @testset "graphs from polytopes" begin
        @test n_vertices(egtriangle) == 3
        @test n_edges(egtriangle) == 3
        @test n_vertices(dgtriangle) == 3
        @test n_edges(dgtriangle) == 3
        @test n_vertices(egcube) == 8
        @test n_edges(egcube) == 12
        @test n_vertices(dgcube) == 6
        @test n_edges(dgcube) == 12

        @test is_isomorphic(dgtriangle, egtriangle)

        @test is_isomorphic(egcr, dgcube)
        @test !is_isomorphic(egcr, egcube)

        # unbounded examples
        @test n_vertices(egpos) == 2
        @test n_edges(egpos) == 1

        @test n_vertices(egpl) == 0
        @test n_edges(egpl) == 0
        @test n_vertices(egplc) == 2
        @test n_edges(egplc) == 1

        @test incidence_matrix(egtriangle) == IncidenceMatrix([[1,2],[1,3],[2,3]])

        @test is_isomorphic(dual_graph(convex_hull([0 0 0; 1 0 0], nothing, [0 1 0])), Graph{Undirected}(2))
        @test is_isomorphic(dual_graph(convex_hull([0 0 0], [0 0 1; 0 1 0; 1 0 0])), complete_graph(3))
        g = dual_graph(convex_hull([0 0 0; 1 0 0], [0 0 1; 0 1 0]))
        @test n_vertices(g) == 4
        @test n_edges(g) == 5
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

    @testset "grap_from_edges" begin
        x1 = [[5,6],[7,8],[11,12]]
        G1 = graph_from_edges(x1)

        @test n_vertices(G1) == 12
        @test n_edges(G1) == 3
      
        x2 = [[11,3],[3,5],[4,5],[2,4],[2,3]]
        G2 = graph_from_edges(Undirected, x2, 13)

        @test n_vertices(G2) == 13
        @test n_edges(G2) == 5

        ei = edges(G2)
        @test length(ei) == 5

        ee = collect(ei)
        @test length(ei) == 0
        @test collect(ei) == Edge[]

        GG2 = graph_from_edges(Undirected, ee, 13)
        @test is_isomorphic(G2, GG2)

    end

    @testset "adjacency_matrix laplacian_matrix" begin
      G0 = Graph{Directed}(3)
      add_edge!(G0,1,2)
      add_edge!(G0,1,3)
      @test matrix(ZZ, adjacency_matrix(G0)) == matrix(ZZ, [0 1 1; 0 0 0; 0 0 0])
      @test laplacian_matrix(G0) == matrix(ZZ, [2 -1 -1; 0 0 0; 0 0 0])
      G1 = vertex_edge_graph(cube(2))
      @test matrix(ZZ, adjacency_matrix(G1)) == matrix(ZZ, [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0])
      @test laplacian_matrix(G1) == matrix(ZZ, [2 -1 -1 0; -1 2 0 -1; -1 0 2 -1; 0 -1 -1 2])
    end
end
