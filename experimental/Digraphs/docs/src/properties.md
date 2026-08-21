```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Properties of Digraphs

The properties in this chapter mirror Chapter 6 of the GAP Digraphs manual. Every property is a boolean-valued function of a digraph.

## Vertex properties

```@docs
has_a_vertex(D::Digraph)
has_no_vertices(D::Digraph)
```

## Edge properties

```@docs
has_loops(D::Digraph)
is_antisymmetric(D::Digraph)
is_bipartite(D::Digraph)
is_complete_bipartite(D::Digraph)
is_complete(D::Digraph)
is_complete_multipartite(D::Digraph)
is_empty(D::Digraph)
is_null(D::Digraph)
is_equivalence(D::Digraph)
is_functional(D::Digraph)
is_permutation(D::Digraph)
is_multi(D::Digraph)
is_nonempty(D::Digraph)
is_reflexive(D::Digraph)
is_symmetric(D::Digraph)
is_tournament(D::Digraph)
is_transitive(D::Digraph)
```

## Edge Weights

```@docs
edge_weights(D::Digraph)
edge_weights_mutable_copy(D::Digraph)
edge_weighted_digraph(D::Digraph, weights::Vector{<:AbstractVector})
edge_weighted_digraph_total_weight(D::Digraph)
edge_weighted_digraph_minimum_spanning_tree(D::Digraph)
edge_weighted_digraph_shortest_paths(D::Digraph)
edge_weighted_digraph_shortest_paths(D::Digraph, source::Int)
edge_weighted_digraph_shortest_path(D::Digraph, s::Int, t::Int)
digraph_maximum_flow(D::Digraph, s::Int, t::Int)
digraph_minimum_cut(D::Digraph, s::Int, t::Int)
digraph_minimum_cut_set(D::Digraph, s::Int, t::Int)
random_unique_edge_weighted_digraph(n::Integer; mut::Bool=false)
is_negative_edge_weighted(D::Digraph)
```

## Orders

```@docs
is_preorder(D::Digraph)
is_quasiorder(D::Digraph)
is_partial_order(D::Digraph)
is_meet_semilattice(D::Digraph)
is_join_semilattice(D::Digraph)
is_lattice(D::Digraph)
digraph_meet_table(D::Digraph)
digraph_join_table(D::Digraph)
is_order_ideal(D::Digraph, subset::Vector{<:Integer})
is_order_filter(D::Digraph, subset::Vector{<:Integer})
is_upper_semimodular(D::Digraph)
is_lower_semimodular(D::Digraph)
is_distributive_lattice(D::Digraph)
is_modular_lattice(D::Digraph)
```

## Regularity

```@docs
is_in_regular(D::Digraph)
is_out_regular(D::Digraph)
is_regular(D::Digraph)
is_distance_regular(D::Digraph)
```

## Connectivity and cycles

```@docs
is_acyclic(D::Digraph)
is_chain(D::Digraph)
is_connected(D::Digraph)
is_biconnected(D::Digraph)
is_bridgeless(D::Digraph)
is_strongly_connected(D::Digraph)
is_aperiodic(D::Digraph)
is_directed_tree(D::Digraph)
is_directed_forest(D::Digraph)
is_undirected_tree(D::Digraph)
is_undirected_forest(D::Digraph)
is_eulerian(D::Digraph)
is_hamiltonian(D::Digraph)
is_cycle(D::Digraph)
```

## Planarity

```@docs
is_planar(D::Digraph)
is_outer_planar(D::Digraph)
```

## Homomorphisms and transformations

```@docs
is_digraph_core(D::Digraph)
is_edge_transitive(D::Digraph)
is_vertex_transitive(D::Digraph)
is_2_edge_transitive(D::Digraph)
is_cayley(D::Digraph)
```
