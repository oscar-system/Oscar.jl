@testset "Graphs" begin
    @testset "core functionality" begin
        g = Graph{Directed}(5)
        @test n_vertices(g) == 5
        @test n_edges(g) == 0
        @test vertices(g) == 1:5
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
        @test vertices(g) == 1:10
        @test rem_vertices!(g, [2, 4, 6, 11])
        @test n_vertices(g) == 7
        @test vertices(g) == 1:7

        @test degree(g) == zeros(Int,7)
        @test degree(g,2) == 0

        g = Graph{Directed}(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        @test signed_incidence_matrix(g) == Matrix([-1 0 1; 1 -1 0; 0 1 -1; 0 0 0])

        @test indegree(g) == [1,1,1,0]
        @test outdegree(g) == [1,1,1,0]

        @test indegree(g, 1) == 1
        @test outdegree(g, 1) == 1

        e = Edge(1,2)
        @test 1 in e
        @test 2 in e
        @test !(3 in e)
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

        @test incidence_matrix(egtriangle) == incidence_matrix([[1,2],[1,3],[2,3]])

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
        @test connectivity(g) == 0
        @test length(connected_components(g)) == 5

        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        add_edge!(g, 4, 5)
        @test !is_connected(g)
        @test connectivity(g) == 0
        @test length(connected_components(g)) == 2

        add_edge!(g, 3, 5)
        @test is_connected(g)
        @test connectivity(g) == 1
        @test length(connected_components(g)) == 1
        @test diameter(g) == 3
    end

    @testset "errors" begin
        g = Graph{Undirected}(1)
        @test !add_edge!(g,1,2)
    end

    @testset "graph_from_edges" begin
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

      @testset "graph_from_labeled_edges" begin
        G1 = graph_from_edges([[1, 2], [3, 4]])
        vertex_labels = Dict(1 => 2, 3 => 4)
        label!(G1, nothing, vertex_labels; name=:color)
        @test_throws ArgumentError G1.color[1, 2]
        @test G1.color[1] == 2

        edge_labels = Dict((1, 2) => 1, (3, 4) => 2)
        label!(G1, edge_labels, nothing; name=:color)
        @test G1.color[1, 2] == 1
        @test G1.color[1] == 2
        
        edge_labels = Dict((5, 6) => 4, (7, 8) => 3)
        G2 = graph_from_labeled_edges(edge_labels)
        @test G2.label[6, 5] == G2.label[5, 6] == 4
        @test_throws ArgumentError G2.label[6, 7]
        @test_throws ArgumentError G2.label[6]
        label!(G2, nothing, Dict(1 => 1))
        @test G2.label[1] ==  1
        @test G2.label[6, 5] == G2.label[5, 6] == 4

        vertex_labels = Dict(9 => 10)
        @test_throws ArgumentError graph_from_labeled_edges(Directed, edge_labels, vertex_labels)

        edge_labels = Dict{NTuple{2, Int}, QQFieldElem}((5, 6) => 3//4, (7, 8) => 3)
        vertex_labels = Dict{Int, QQFieldElem}(3 => 1//4, 9 => 10)
        G3 = graph_from_labeled_edges(Directed, edge_labels, vertex_labels; n_vertices=9)
        @test_throws ArgumentError G3.label[10]
        @test_throws ArgumentError G3.label[6, 5]
        @test G3.label[7, 8] == 3
        @test G3.label[9] == 10
        @test G3.label[1] == 0
        @test labelings(G3) == [:label]
        @test G3.label[5, 6] isa QQFieldElem
        @test G3.label[3] isa QQFieldElem
        
        vertex_labels[5] = 3
        edge_labels[2, 3] = 4
        @test_throws ArgumentError label!(G1, nothing, vertex_labels)
        @test_throws ArgumentError label!(G1, edge_labels, vertex_labels)
        @test_throws ArgumentError label!(G1, edge_labels, nothing)
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

    @testset "is_bipartite" begin
      G0 = Graph{Undirected}(3)
      add_edge!(G0,1,2)
      add_edge!(G0,1,3)
      add_edge!(G0,2,3)
      @test is_bipartite(G0) == false

      G1 = graph_from_edges([[1,2],[2,3],[3,4]])
      @test is_bipartite(G1) == true
    end

    @testset "maximal_cliques" begin
      G = complete_bipartite_graph(2, 2)
      @test maximal_cliques(G) == Set{Set{Int}}(Set.([[1, 3], [1, 4], [2, 3], [2, 4]]))
    end

    @testset "is_acyclic" begin
      G = graph_from_edges(Directed, [[1, 2], [2, 3], [3, 1]])
      @test !is_acyclic(G)
      rem_edge!(G, 3, 1)
      @test is_acyclic(G)
    end

    @testset "subgraph" begin
      G = graph_from_edges(Directed, [[1, 2], [2, 3], [3, 1]])
      sg = induced_subgraph(G, [1, 2])
      @test ne(sg) == 1
      @test nv(sg) == 2
      @test sg isa Graph{Directed}
      G2 = complete_bipartite_graph(3, 3)
      sg2 = induced_subgraph(G2, [1, 2, 3])
      @test ne(sg2) == 0
      @test nv(sg2) == 3

      G3 = graph_from_labeled_edges(Undirected, Dict((1, 2) => 4, (2, 3) => 5, (1, 3) => 6), Dict(3 => 9); name=:color)
      sg3 = induced_subgraph(G3, [3, 2])
      @test ne(sg3) == 1
      @test nv(sg3) == 2
      @test sg3.color[2] == 9
      @test sg3.vertexlabels[2] == 3
      @test sg3.color[1,2] == G3.color[2,3] == 5

      G4 = complete_graph(4)
      label!(G4,nothing,Dict(1=>"first",2=>"second",3=>"third",4=>"fourth"), name=:vertexlabels)
      sg4 = induced_subgraph(G4, [2,4])
      @test sg4.vertexlabels[1] == "second"
      @test sg4.vertexlabels[2] == "fourth"
    end
end
