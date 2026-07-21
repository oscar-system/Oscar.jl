@testset "Digraphs" begin

  @testset "Constructors" begin
    D = cycle_digraph(4)
    @test nv(D) == 4
    @test ne(D) == 4

    D2 = complete_digraph(3)
    @test nv(D2) == 3
    @test ne(D2) == 6

    D3 = null_digraph()
    @test nv(D3) == 0
    @test ne(D3) == 0

    D4 = complete_bipartite_digraph(2, 3)
    @test nv(D4) == 5

    D5 = digraph_from_edges(4, [(1,2), (2,3), (3,4), (4,1)])
    @test nv(D5) == 4
    @test ne(D5) == 4

    D6 = digraph_from_adjacency_matrix([0 1 0; 0 0 1; 1 0 0])
    @test nv(D6) == 3
    @test ne(D6) == 3
  end

  @testset "Properties" begin
    D = cycle_digraph(4)
    @test is_connected(D)
    @test is_strongly_connected(D)
    @test !is_acyclic(D)
    @test is_bipartite(D)
    @test !is_complete(D)
    @test is_regular(D)

    K4 = complete_digraph(4)
    @test is_complete(K4)
    @test is_symmetric(K4)
    @test !is_antisymmetric(K4)
    @test is_vertex_transitive(K4)
    @test is_edge_transitive(K4)

    D3 = digraph_from_edges(3, [(1,2), (2,3)])
    @test is_acyclic(D3)
    @test is_directed_tree(D3)
    @test !is_strongly_connected(D3)
  end

  @testset "Attributes" begin
    D = cycle_digraph(4)
    @test out_neighbours(D) == [[2], [3], [4], [1]]
    @test in_neighbours(D) == [[4], [1], [2], [3]]
    @test out_degrees(D) == [1, 1, 1, 1]
    @test in_degrees(D) == [1, 1, 1, 1]
    @test has_edge(D, 1, 2)
    @test !has_edge(D, 2, 1)
    @test has_vertex(D, 1)
    @test !has_vertex(D, 5)

    K3 = complete_digraph(3)
    @test chromatic_number(K3) == 3
    @test clique_number(K3) == 3
  end

  @testset "Operations" begin
    R = reverse_digraph(cycle_digraph(4))
    @test out_neighbours(R) == [[4], [1], [2], [3]]

    U = digraph_disjoint_union(cycle_digraph(3), cycle_digraph(4))
    @test nv(U) == 7

    CP = digraph_cartesian_product(cycle_digraph(3), cycle_digraph(3))
    @test nv(CP) == 9
  end

  @testset "Isomorphisms" begin
    D1 = cycle_digraph(4)
    D2 = digraph_from_edges(4, [(1,2), (2,3), (3,4), (4,1)])
    @test is_isomorphic(D1, D2)

    K4 = complete_digraph(4)
    aut = automorphism_group(K4)
    @test GAPWrap.Size(aut) == 24
  end

end
