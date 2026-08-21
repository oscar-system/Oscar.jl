```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Attributes of Digraphs

The attributes and operations in this chapter mirror Chapter 5 of the GAP Digraphs manual. A few functions from other chapters of the manual (cliques and input/output) are collected here as well, because they live in the same source file.

## Vertices and edges

```@docs
nv(D::Digraph)
ne(D::Digraph)
digraph_vertices(D::Digraph)
digraph_edges(D::Digraph)
digraph_nr_adjacencies(D::Digraph)
digraph_nr_adjacencies_without_loops(D::Digraph)
digraph_nr_loops(D::Digraph)
digraph_sinks(D::Digraph)
digraph_sources(D::Digraph)
digraph_topological_sort(D::Digraph)
digraph_vertex_label(D::Digraph, i::Int)
set_digraph_vertex_label!(D::Digraph, i::Int, obj)
digraph_vertex_labels(D::Digraph)
set_digraph_vertex_labels!(D::Digraph, list)
digraph_edge_label(D::Digraph, i::Int, j::Int)
set_digraph_edge_label!(D::Digraph, i::Int, j::Int, obj)
digraph_edge_labels(D::Digraph)
set_digraph_edge_labels!(D::Digraph, labels)
digraph_in_edges(D::Digraph, v::Int)
digraph_out_edges(D::Digraph, v::Int)
is_digraph_edge(D::Digraph, e::Vector{<:Integer})
is_matching(D::Digraph, edges::Vector)
is_maximal_matching(D::Digraph, edges::Vector)
is_maximum_matching(D::Digraph, edges::Vector)
is_perfect_matching(D::Digraph, edges::Vector)
digraph_maximal_matching(D::Digraph)
digraph_maximum_matching(D::Digraph)
has_edge(D::Digraph, s::Int, t::Int)
has_vertex(D::Digraph, v::Int)
```

## Neighbours and degree

```@docs
adjacency_matrix(D::Digraph)
adjacency_matrix_mutable_copy(D::Digraph)
characteristic_polynomial(D::Digraph)
boolean_adjacency_matrix(D::Digraph)
boolean_adjacency_matrix_mutable_copy(D::Digraph)
digraph_adjacency_function(D::Digraph)
digraph_source(D::Digraph)
digraph_range(D::Digraph)
out_neighbours(D::Digraph)
out_neighbors(D::Digraph)
out_neighbours_mutable_copy(D::Digraph)
out_neighbors_mutable_copy(D::Digraph)
in_neighbours(D::Digraph)
in_neighbors(D::Digraph)
in_neighbours_mutable_copy(D::Digraph)
in_neighbors_mutable_copy(D::Digraph)
out_degrees(D::Digraph)
out_degree_sequence(D::Digraph)
out_degree_set(D::Digraph)
in_degrees(D::Digraph)
in_degree_sequence(D::Digraph)
in_degree_set(D::Digraph)
out_degree_of_vertex(D::Digraph, v::Int)
out_neighbours_of_vertex(D::Digraph, v::Int)
out_neighbors_of_vertex(D::Digraph, v::Int)
in_degree_of_vertex(D::Digraph, v::Int)
in_neighbours_of_vertex(D::Digraph, v::Int)
in_neighbors_of_vertex(D::Digraph, v::Int)
digraph_loops(D::Digraph)
degree_matrix(D::Digraph)
laplacian_matrix(D::Digraph)
```

## Orders

```@docs
partial_order_digraph_meet_of_vertices(D::Digraph, u::Int, v::Int)
partial_order_digraph_join_of_vertices(D::Digraph, u::Int, v::Int)
non_upper_semimodular_pair(D::Digraph)
non_lower_semimodular_pair(D::Digraph)
```

## Reachability and connectivity

```@docs
digraph_diameter(D::Digraph)
digraph_shortest_distance(D::Digraph, s::Int, t::Int)
digraph_shortest_distance(D::Digraph, list::Vector{<:Integer})
digraph_shortest_distance(D::Digraph, list1::Vector{<:Integer}, list2::Vector{<:Integer})
shortest_distances(D::Digraph)
digraph_longest_distance_from_vertex(D::Digraph, v::Int)
digraph_distance_set(D::Digraph, v::Int, distance::Int)
digraph_girth(D::Digraph)
digraph_odd_girth(D::Digraph)
digraph_undirected_girth(D::Digraph)
digraph_connected_components(D::Digraph)
digraph_nr_connected_components(D::Digraph)
digraph_connected_component(D::Digraph, v::Int)
digraph_strongly_connected_components(D::Digraph)
digraph_nr_strongly_connected_components(D::Digraph)
digraph_strongly_connected_component(D::Digraph, v::Int)
digraph_bicomponents(D::Digraph)
articulation_points(D::Digraph)
minimal_cyclic_edge_cut(D::Digraph)
bridges(D::Digraph)
strong_orientation(D::Digraph)
strong_orientation_attr(D::Digraph)
digraph_period(D::Digraph)
digraph_floyd_warshall(D::Digraph, func, nopath, edge)
is_reachable(D::Digraph, s::Int, t::Int)
is_digraph_path(D::Digraph, v::Vector{<:Integer}, a::Vector{<:Integer})
vertices_reachable_from(D::Digraph, root::Int)
digraph_path(D::Digraph, s::Int, t::Int)
digraph_shortest_path(D::Digraph, s::Int, t::Int)
digraph_random_walk(D::Digraph, v::Int, t::Int)
digraph_absorption_probabilities(D::Digraph)
digraph_absorption_expected_steps(D::Digraph)
dominators(D::Digraph, root::Int)
dominator_tree(D::Digraph, root::Int)
iterator_of_paths(D::Digraph, u::Int, v::Int)
digraph_all_simple_circuits(D::Digraph)
digraph_longest_simple_circuit(D::Digraph)
digraph_all_undirected_simple_circuits(D::Digraph)
digraph_all_chordless_cycles(D::Digraph)
digraph_all_chordless_cycles_of_maximal_length(D::Digraph, max_length::Int)
facial_walks(D::Digraph, list::Vector{<:Vector{<:Integer}})
digraph_layers(D::Digraph, v::Int)
digraph_degeneracy(D::Digraph)
digraph_degeneracy_ordering(D::Digraph)
hamiltonian_path(D::Digraph)
nr_spanning_trees(D::Digraph)
digraph_dijkstra(D::Digraph, s::Int)
digraph_vertex_connectivity(D::Digraph)
digraph_cycle_basis(D::Digraph)
digraph_is_king(D::Digraph, v::Int, k::Int)
digraph_kings(D::Digraph, n::Int)
```

## Cayley graphs of groups

```@docs
group_of_cayley_digraph(D::Digraph)
semigroup_of_cayley_digraph(D::Digraph)
generators_of_cayley_digraph(D::Digraph)
```

## Associated semigroups

```@docs
as_semigroup(filt, D::Digraph)
as_monoid(filt, D::Digraph)
as_semigroup(filt, Y::Digraph, gps, homs)
```

## Planarity

```@docs
kuratowski_planar_subdigraph(D::Digraph)
kuratowski_outer_planar_subdigraph(D::Digraph)
planar_embedding(D::Digraph)
outer_planar_embedding(D::Digraph)
subdigraph_homeomorphic_to_k23(D::Digraph)
subdigraph_homeomorphic_to_k33(D::Digraph)
subdigraph_homeomorphic_to_k4(D::Digraph)
dual_planar_graph(D::Digraph)
```

## Hashing

```@docs
digraph_hash(D::Digraph)
```

## Finding cliques

```@docs
clique_number(D::Digraph)
```

## Reading and writing digraphs to a file

```@docs
graph6_string(D::Digraph)
sparse6_string(D::Digraph)
```
