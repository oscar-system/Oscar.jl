@testset "Digraphs" begin

  @testset "3.1 Creating digraphs" begin
    # adjacency list, labels, integer source/range, function, named digraph
    d = digraph([[2, 5, 8, 10], [2, 3, 4, 2, 5, 6, 8, 9, 10], [1],
                 [3, 5, 7, 8, 10], [2, 5, 7], [3, 6, 7, 9, 10], [1, 4],
                 [1, 5, 9], [1, 2, 7, 8], [3, 5]])
    @test nv(d) == 10
    @test ne(d) == 38
    @test is_digraph(d)
    @test is_immutable_digraph(d)
    @test !is_mutable_digraph(d)

    d = digraph(["a", "b", "c"], ["a"], ["b"])
    @test nv(d) == 3 && ne(d) == 1

    d = digraph(5, [1, 2, 2, 4, 1, 1], [2, 3, 5, 5, 1, 1])
    @test nv(d) == 5 && ne(d) == 6

    d = digraph(collect(1:10), (x, y) -> true)
    @test nv(d) == 10 && ne(d) == 100

    d = digraph("Diamond")
    @test nv(d) == 4 && ne(d) == 10

    d = digraph([[2], [1]]; mut=true)
    @test is_mutable_digraph(d)

    # adjacency matrix
    d = digraph_from_adjacency_matrix([0 1 0; 1 0 1; 0 1 0])
    @test nv(d) == 3 && ne(d) == 4
    d = digraph_from_adjacency_matrix([true false; false true])
    @test nv(d) == 2 && ne(d) == 2

    # edges
    d = digraph_from_edges([[1, 2], [2, 3], [3, 1]])
    @test nv(d) == 3 && ne(d) == 3
    d = digraph_from_edges([(1, 2), (2, 3), (3, 1)], 5)
    @test nv(d) == 5 && ne(d) == 3

    # edge orbits
    d = edge_orbit_digraph(symmetric_group(4), [[1, 2], [3, 4]])
    @test nv(d) == 4 && ne(d) == 12
    d = edge_orbit_digraph(symmetric_group(4), [(1, 2), (3, 4)])
    @test nv(d) == 4 && ne(d) == 12

    # in-neighbours
    d = digraph_from_in_neighbours([[2], [1, 3], [2]])
    @test nv(d) == 3 && ne(d) == 4
    d2 = digraph_from_in_neighbors([[2], [1, 3], [2]])
    @test nv(d2) == 3 && ne(d2) == 4

    # Cayley digraph
    d = cayley_digraph(symmetric_group(4))
    @test nv(d) == 24 && ne(d) == 48
    @test is_cayley_digraph(d)

    # named digraph database
    @test length(list_named_digraphs("diamond")) == 1
    @test length(list_named_digraphs("dia")) >= 1
  end

  @testset "3.2 Changing representations" begin
    D = digraph([[3, 2], [1, 2], [2], [3, 4]])
    br = as_binary_relation(D)
    @test isa(br, GapObj)
    D2 = as_digraph(br)
    @test nv(D2) == 4 && ne(D2) == 7

    f = GAP.Globals.Transformation(GapObj([4, 3, 3, 1, 7, 9, 10, 4, 2, 3]))
    D = as_digraph(f)
    @test nv(D) == 10 && ne(D) == 10
    @test_throws ArgumentError as_digraph(f, 5)

    D = digraph([[1], [3], [2]])
    t = as_transformation(D)
    @test isa(t, GapObj)
    @test_throws ArgumentError as_transformation(digraph([[2, 3], [], []]))

    @test isa(to_grape_graph(cycle_digraph(4)), GapObj)
    @test isa(as_grape_graph(cycle_digraph(4)), GapObj)
  end

  @testset "3.3 New digraphs from old" begin
    # copies
    d = cycle_digraph(4)
    @test nv(digraph_immutable_copy(d)) == 4
    @test is_mutable_digraph(digraph_mutable_copy(d))
    @test is_immutable_digraph(digraph_copy_same_mutability(d))
    @test is_immutable_digraph(digraph_copy(d))
    @test is_immutable_digraph(digraph_immutable_copy_if_immutable(d))
    @test is_immutable_digraph(digraph_immutable_copy_if_mutable(digraph_mutable_copy(d)))
    m = digraph_mutable_copy(d)
    @test is_mutable_digraph(digraph_mutable_copy_if_mutable(m))
    @test is_mutable_digraph(digraph_mutable_copy_if_immutable(d))

    # induced and reduced subdigraphs
    d = digraph([[1, 1, 2, 3, 4, 4], [1, 3, 4], [3, 1], [1, 1]])
    s = induced_subdigraph(d, [1, 3, 4])
    @test nv(s) == 3 && ne(s) == 9
    d = digraph([[1, 2], [], [], [1, 4], []])
    @test nv(reduced_digraph(d)) == 3 && ne(reduced_digraph(d)) == 4
    @test nv(reduced_digraph_attr(d)) == 3

    # symmetric / antisymmetric subdigraphs
    d = digraph([[2, 2], [1, 3], [4], [3, 1]])
    @test nv(maximal_symmetric_subdigraph(d)) == 4 && ne(maximal_symmetric_subdigraph(d)) == 4
    @test ne(maximal_symmetric_subdigraph_attr(d)) == 4
    @test ne(maximal_anti_symmetric_subdigraph(d)) == 4
    @test ne(maximal_anti_symmetric_subdigraph_attr(d)) == 4
    d = digraph([[1, 1, 2], [1], [1, 2]])
    @test ne(maximal_symmetric_subdigraph_without_loops(d)) == 2
    @test ne(maximal_symmetric_subdigraph_without_loops_attr(d)) == 2

    # spanning forests and trees
    d = digraph([[1, 2, 1, 3], [1], [4], [3, 4, 3]])
    @test nv(undirected_spanning_forest(d)) == 4
    @test nv(undirected_spanning_forest_attr(d)) == 4
    d = complete_digraph(4)
    @test nv(undirected_spanning_tree(d)) == 4 && ne(undirected_spanning_tree(d)) == 6
    @test nv(undirected_spanning_tree_attr(d)) == 4
    @test_throws ArgumentError undirected_spanning_forest(empty_digraph(0))
    @test_throws ArgumentError undirected_spanning_tree(digraph([[2], [1], []]))

    # shortest path spanning tree
    @test nv(digraph_shortest_path_spanning_tree(complete_digraph(4), 1)) == 4
    @test_throws ArgumentError digraph_shortest_path_spanning_tree(digraph([[2], [3], []]), 3)

    # quotient digraph
    d = digraph([[2, 1], [4], [1], [1, 3, 4]])
    q = quotient_digraph(d, [[1], [2, 4], [3]])
    @test nv(q) == 3 && ne(q) == 6

    # reverse, dual, closures, reductions
    d = digraph([[3], [1, 3, 5], [1], [1, 2, 4], [2, 3, 5]])
    @test ne(digraph_reverse(d)) == 11
    @test ne(digraph_reverse_attr(d)) == 11
    d = digraph([[2, 3], [], [4, 6], [5], [], [7, 8, 9], [], [], []])
    @test nv(digraph_dual(d)) == 9 && ne(digraph_dual(d)) == 73
    @test ne(digraph_dual_attr(d)) == 73
    d = digraph([[1, 2, 3], [2, 4], [1], [3, 4]])
    @test ne(digraph_symmetric_closure(d)) == 11
    @test ne(digraph_symmetric_closure_attr(d)) == 11
    d = digraph([[2], [3], []])
    @test ne(digraph_transitive_closure(d)) == 3
    @test ne(digraph_transitive_closure_attr(d)) == 3
    @test ne(digraph_reflexive_transitive_closure(d)) == 6
    @test ne(digraph_reflexive_transitive_closure_attr(d)) == 6
    d = digraph([[2, 3], [3], []])
    @test ne(digraph_transitive_reduction(d)) == 2
    @test ne(digraph_transitive_reduction_attr(d)) == 2
    @test ne(digraph_reflexive_transitive_reduction(d)) == 2
    @test ne(digraph_reflexive_transitive_reduction_attr(d)) == 2

    # adding vertices and edges
    d = complete_digraph(3)
    @test nv(digraph_add_vertex(d)) == 4
    @test nv(digraph_add_vertices(d, 3)) == 6
    @test nv(digraph_add_vertices(empty_digraph(2), ["a", "b"])) == 4
    d = digraph([[2], [3], []])
    @test ne(digraph_add_edge(d, [3, 1])) == 3
    @test ne(digraph_add_edge(d, 1, 2)) == 3
    @test ne(digraph_add_edges(empty_digraph(3), [[1, 2], [2, 3]])) == 2
    @test ne(digraph_add_edges(empty_digraph(3), [(1, 2), (2, 3)])) == 2
    @test ne(digraph_add_all_loops(empty_digraph(5))) == 5
    @test ne(digraph_add_all_loops_attr(empty_digraph(3))) == 3

    # removing vertices and edges
    d = digraph(["a", "b", "c"], ["a", "a", "b", "c", "c"], ["b", "c", "a", "a", "c"])
    @test nv(digraph_remove_vertex(d, 2)) == 2 && ne(digraph_remove_vertex(d, 2)) == 3
    d = digraph([[3], [1, 3, 5], [1], [1, 2, 4], [2, 3, 5]])
    @test nv(digraph_remove_vertices(d, [2, 4])) == 3 && ne(digraph_remove_vertices(d, [2, 4])) == 4
    d = cycle_digraph(4)
    @test ne(digraph_remove_edge(d, [4, 1])) == 3
    @test ne(digraph_remove_edge(complete_digraph(2), 1, 2)) == 1
    @test ne(digraph_remove_edges(complete_digraph(3), [[1, 2], [2, 1]])) == 4
    d = digraph([[1, 2, 4], [1, 4], [3, 4], [1, 4, 5], [1, 5]])
    @test ne(digraph_remove_loops(d)) == 8
    @test ne(digraph_remove_loops_attr(d)) == 8
    d = digraph([[1, 2, 3, 2], [1, 1, 3], [2, 2, 2]])
    @test ne(digraph_remove_all_multiple_edges(d)) == 6
    @test ne(digraph_remove_all_multiple_edges_attr(d)) == 6

    # contracting and reversing edges
    d = digraph_from_edges([[1, 2], [2, 1]])
    @test nv(digraph_contract_edge(d, 1, 2)) == 1 && ne(digraph_contract_edge(d, 1, 2)) == 0
    d = digraph([[2], [3], []])
    @test ne(digraph_reverse_edges(d, [[1, 2]])) == 2
    @test ne(digraph_reverse_edge(d, [1, 2])) == 2
    @test ne(digraph_reverse_edge(d, 2, 3)) == 2

    # unions, joins, products
    @test nv(digraph_disjoint_union(cycle_digraph(3), cycle_digraph(4))) == 7
    @test nv(digraph_disjoint_union([cycle_digraph(3), cycle_digraph(4)])) == 7
    @test nv(digraph_disjoint_union(cycle_digraph(3), cycle_digraph(4), cycle_digraph(5))) == 12
    d1 = digraph([[2], [1]])
    d2 = digraph([[2, 3], [2], [1]])
    @test nv(digraph_edge_union(d1, d2)) == 3 && ne(digraph_edge_union(d1, d2)) == 6
    @test ne(digraph_join(complete_digraph(3), cycle_digraph(3))) == 27
    @test nv(digraph_cartesian_product(chain_digraph(4), cycle_digraph(3))) == 12 &&
          ne(digraph_cartesian_product(chain_digraph(4), cycle_digraph(3))) == 21
    @test nv(digraph_cartesian_product([chain_digraph(4), cycle_digraph(3)])) == 12
    @test nv(digraph_direct_product(chain_digraph(4), cycle_digraph(3))) == 12 &&
          ne(digraph_direct_product(chain_digraph(4), cycle_digraph(3))) == 9

    s3 = digraph_symmetric_closure(cycle_digraph(3))
    s4 = digraph_symmetric_closure(cycle_digraph(4))
    @test nv(conormal_product(s3, s4)) == 12 && ne(conormal_product(s3, s4)) == 120
    @test nv(homomorphic_product(petersen_digraph(),
                                 digraph_symmetric_closure(chain_digraph(4)))) == 40 &&
          ne(homomorphic_product(petersen_digraph(),
                                 digraph_symmetric_closure(chain_digraph(4)))) == 460
    @test nv(lexicographic_product(s3, digraph_symmetric_closure(chain_digraph(2)))) == 6 &&
          ne(lexicographic_product(s3, digraph_symmetric_closure(chain_digraph(2)))) == 30
    @test nv(modular_product(digraph([[1], [1, 2]]), digraph([Int[], [2]]))) == 4 &&
          ne(modular_product(digraph([[1], [1, 2]]), digraph([Int[], [2]]))) == 4
    @test nv(strong_product(digraph_symmetric_closure(chain_digraph(3)),
                            digraph_symmetric_closure(chain_digraph(4)))) == 12 &&
          ne(strong_product(digraph_symmetric_closure(chain_digraph(3)),
                            digraph_symmetric_closure(chain_digraph(4)))) == 58

    # product aliases
    @test nv(digraph_lexicographic_product(s3, s4)) == nv(lexicographic_product(s3, s4))
    @test nv(digraph_modular_product(digraph([[1], [1, 2]]), digraph([Int[], [2]]))) == 4
    @test nv(digraph_strong_product(chain_digraph(3), chain_digraph(4))) ==
          nv(strong_product(chain_digraph(3), chain_digraph(4)))

    # projections
    cp = digraph_cartesian_product(chain_digraph(3), cycle_digraph(4))
    @test length(digraph_cartesian_product_projections(cp)) == 2
    dp = digraph_direct_product(chain_digraph(3), cycle_digraph(4))
    @test length(digraph_direct_product_projections(dp)) == 2

    # line digraphs and double digraphs
    @test nv(line_digraph(complete_digraph(3))) == 6 && ne(line_digraph(complete_digraph(3))) == 12
    @test nv(edge_digraph(complete_digraph(3))) == 6
    @test nv(line_undirected_digraph(complete_digraph(3))) == 3 &&
          ne(line_undirected_digraph(complete_digraph(3))) == 6
    @test nv(edge_undirected_digraph(complete_digraph(3))) == 3
    @test nv(double_digraph(digraph([[2], [3], [1]]))) == 6 &&
          ne(double_digraph(digraph([[2], [3], [1]]))) == 12
    @test nv(bipartite_double_digraph(digraph([[2], [3], [1]]))) == 6 &&
          ne(bipartite_double_digraph(digraph([[2], [3], [1]]))) == 6

    # distance digraph, closure, Mycielskian
    d = cycle_digraph(5)
    @test ne(distance_digraph(d, 2)) == 5
    @test ne(distance_digraph(d, [1, 2])) == 10
    d = digraph_remove_edges(complete_digraph(6), [[1, 2], [2, 1]])
    @test ne(digraph_closure(d, 6)) == 30
    @test nv(digraph_mycielskian(cycle_digraph(2))) == 5 &&
          ne(digraph_mycielskian(cycle_digraph(2))) == 10
    @test ne(digraph_mycielskian_attr(cycle_digraph(2))) == 10

    # in-place semantics for mutable inputs, new objects for immutable inputs
    m = digraph([[2], [3], []]; mut=true)
    m2 = digraph_reverse(m)
    @test is_mutable_digraph(m2)
    @test GapObj(m2) === GapObj(m)
    @test out_neighbours(m) == [[], [1], [2]]

    im = digraph([[2], [3], []])
    im2 = digraph_reverse(im)
    @test is_immutable_digraph(im2)
    @test GapObj(im2) !== GapObj(im)
    @test out_neighbours(im) == [[2], [3], []]

    m = digraph([[2], [1]]; mut=true)
    m2 = digraph_add_vertex(m)
    @test nv(m2) == 3
    @test is_mutable_digraph(m2)
    @test GapObj(m2) === GapObj(m)
  end

  @testset "3.4 Random digraphs" begin
    d = random_digraph(20)
    @test nv(d) == 20
    @test is_immutable_digraph(d)
    d = random_digraph(20, 0.5)
    @test nv(d) == 20
    d = random_digraph(20; mut=:mut)
    @test is_mutable_digraph(d)
    d = random_digraph(20; mut=:connected)
    @test is_connected(d)
    d = random_digraph(20, 0.5; mut=:symmetric)
    @test is_symmetric(d)
    d = random_multi_digraph(20, 30)
    @test nv(d) == 20 && ne(d) == 30
    d = random_tournament(20)
    @test nv(d) == 20 && is_tournament(d)
    d = random_lattice(10)
    @test nv(d) >= 10 && nv(d) <= 20
    @test is_lattice(d)
  end

  @testset "3.5 Standard examples" begin
    @test nv(andrasfai_graph(4)) == 11 && ne(andrasfai_graph(4)) == 44
    @test nv(banana_tree(2, 4)) == 9 && ne(banana_tree(2, 4)) == 16
    @test nv(binary_tree(4)) == 15 && ne(binary_tree(4)) == 14
    @test nv(binomial_tree_graph(4)) == 4 && ne(binomial_tree_graph(4)) == 6
    @test nv(bishops_graph(4, 4)) == 16 && ne(bishops_graph(4, 4)) == 56
    @test nv(bishops_graph("dark", 3, 5)) == 8 && ne(bishops_graph("dark", 3, 5)) == 24
    @test nv(bishop_graph(3, 3)) == nv(bishops_graph(3, 3))
    @test nv(bondy_graph(1)) == 22 && ne(bondy_graph(1)) == 66
    @test nv(book_graph(2)) == 6 && ne(book_graph(2)) == 14
    @test nv(burnt_pancake_graph(3)) == 48 && ne(burnt_pancake_graph(3)) == 144
    @test nv(pancake_graph(4)) == 24 && ne(pancake_graph(4)) == 72
    @test nv(stacked_book_graph(1, 2)) == 4 && ne(stacked_book_graph(1, 2)) == 8
    @test nv(chain_digraph(4)) == 4 && ne(chain_digraph(4)) == 3
    @test nv(circulant_graph(6, [2, 3])) == 6 && ne(circulant_graph(6, [2, 3])) == 18
    @test nv(complete_digraph(3)) == 3 && ne(complete_digraph(3)) == 6
    @test nv(complete_bipartite_digraph(2, 3)) == 5 && ne(complete_bipartite_digraph(2, 3)) == 12
    @test nv(complete_multipartite_digraph([2, 3])) == 5 &&
          ne(complete_multipartite_digraph([2, 3])) == 12
    @test nv(cycle_digraph(4)) == 4 && ne(cycle_digraph(4)) == 4
    @test nv(digraph_cycle(4)) == 4
    @test nv(cycle_symmetric_digraph(7)) == 7 && ne(cycle_symmetric_digraph(7)) == 14
    @test nv(empty_digraph(5)) == 5 && ne(empty_digraph(5)) == 0
    @test nv(null_digraph(5)) == 5
    @test nv(gear_graph(4)) == 9 && ne(gear_graph(4)) == 24
    @test nv(haar_graph(4)) == 6 && ne(haar_graph(4)) == 6
    @test nv(halved_cube_graph(3)) == 4 && ne(halved_cube_graph(3)) == 12
    @test nv(hanoi_graph(4)) == 81 && ne(hanoi_graph(4)) == 240
    @test nv(helm_graph(4)) == 9 && ne(helm_graph(4)) == 24
    @test nv(hypercube_graph(4)) == 16 && ne(hypercube_graph(4)) == 64
    @test nv(johnson_digraph(4, 2)) == 6 && ne(johnson_digraph(4, 2)) == 24
    @test nv(keller_graph(3)) == 64 && ne(keller_graph(3)) == 2176
    @test nv(kings_graph(4, 4)) == 16 && ne(kings_graph(4, 4)) == 84
    @test nv(kneser_graph(5, 2)) == 10 && ne(kneser_graph(5, 2)) == 30
    @test nv(knights_graph(4, 4)) == 16 && ne(knights_graph(4, 4)) == 48
    @test nv(lindgren_sousselier_graph(4)) == 28 && ne(lindgren_sousselier_graph(4)) == 90
    @test nv(lollipop_graph(4, 4)) == 8 && ne(lollipop_graph(4, 4)) == 20
    @test nv(mobius_ladder_graph(4)) == 8 && ne(mobius_ladder_graph(4)) == 24
    @test nv(mycielski_graph(4)) == 11 && ne(mycielski_graph(4)) == 40
    @test nv(odd_graph(4)) == 35 && ne(odd_graph(4)) == 140
    @test nv(path_graph(4)) == 4 && ne(path_graph(4)) == 6
    @test nv(permutation_star_graph(4, 3)) == 24 && ne(permutation_star_graph(4, 3)) == 72
    @test nv(petersen_digraph()) == 10 && ne(petersen_digraph()) == 30
    @test nv(generalised_petersen_graph(7, 2)) == 14 && ne(generalised_petersen_graph(7, 2)) == 42
    @test nv(prism_graph(4)) == 8 && ne(prism_graph(4)) == 24
    @test nv(stacked_prism_graph(5, 2)) == 10 && ne(stacked_prism_graph(5, 2)) == 30
    @test nv(queens_graph(4, 4)) == 16 && ne(queens_graph(4, 4)) == 152
    @test nv(queen_graph(3, 3)) == nv(queens_graph(3, 3))
    @test nv(rooks_graph(4, 4)) == 16 && ne(rooks_graph(4, 4)) == 96
    @test nv(rook_graph(3, 3)) == nv(rooks_graph(3, 3))
    @test nv(square_grid_graph(4, 4)) == 16 && ne(square_grid_graph(4, 4)) == 48
    @test nv(grid_graph(3, 3)) == nv(square_grid_graph(3, 3))
    @test nv(triangular_grid_graph(3, 3)) == 9 && ne(triangular_grid_graph(3, 3)) == 32
    @test nv(star_graph(5)) == 5 && ne(star_graph(5)) == 8
    @test nv(tadpole_graph(4, 4)) == 8 && ne(tadpole_graph(4, 4)) == 16
    @test nv(walsh_hadamard_graph(4)) == 32 && ne(walsh_hadamard_graph(4)) == 256
    @test nv(web_graph(5)) == 15 && ne(web_graph(5)) == 40
    @test nv(wheel_graph(8)) == 8 && ne(wheel_graph(8)) == 28
    @test nv(windmill_graph(4, 3)) == 10 && ne(windmill_graph(4, 3)) == 36

    @test is_mutable_digraph(complete_digraph(3; mut=true))
    @test is_immutable_digraph(complete_digraph(3))
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

    D3 = digraph(3, [1, 2], [2, 3])
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

  @testset "Isomorphisms" begin
    D1 = cycle_digraph(4)
    D2 = digraph(4, [1, 2, 3, 4], [2, 3, 4, 1])
    @test is_isomorphic(D1, D2)

    K4 = complete_digraph(4)
    aut = automorphism_group(K4)
    @test order(aut) == 24
  end

  @testset "Edge-weighted digraph (appendix)" begin
    d = edge_weighted_digraph(digraph([[2], [1]]), [[5], [10]])
    @test nv(d) == 2 && ne(d) == 2
  end

end
