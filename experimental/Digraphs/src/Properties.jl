@doc raw"""
    has_loops(D::Digraph) -> Bool

Return `true` if the digraph `D` has any loops (edges from a vertex to itself).

# Examples
```jldoctest
julia> d = digraph([[1, 2], [2]])
Digraph with 2 vertices, 3 edges

julia> has_loops(d)
true

julia> d2 = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> has_loops(d2)
false
```
"""
has_loops(D::Digraph) = DigraphWrap.DigraphHasLoops(GapObj(D))::Bool

@doc raw"""
    has_no_vertices(D::Digraph) -> Bool

Return `true` if the digraph `D` has no vertices.

# Examples
```jldoctest
julia> d = empty_digraph(0)
Digraph with 0 vertices, 0 edges

julia> has_no_vertices(d)
true
```
"""
has_no_vertices(D::Digraph) = DigraphWrap.DigraphHasNoVertices(GapObj(D))::Bool

@doc raw"""
    is_connected(D::Digraph) -> Bool

Return `true` if the digraph `D` is connected (the underlying undirected
graph is connected).

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges

julia> is_connected(d)
true

julia> d2 = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_connected(d2)
true

```
"""
is_connected(D::Digraph) = DigraphWrap.IsConnectedDigraph(GapObj(D))::Bool

@doc raw"""
    is_strongly_connected(D::Digraph) -> Bool

Return `true` if the digraph `D` is strongly connected (there is a directed
path between every ordered pair of distinct vertices).

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> is_strongly_connected(d)
true

julia> d2 = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_strongly_connected(d2)
false
```
"""
is_strongly_connected(D::Digraph) = DigraphWrap.IsStronglyConnectedDigraph(GapObj(D))::Bool

@doc raw"""
    is_acyclic(D::Digraph) -> Bool

Return `true` if the digraph `D` is acyclic (has no directed cycles).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_acyclic(d)
true

julia> d2 = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> is_acyclic(d2)
false
```
"""
is_acyclic(D::Digraph) = DigraphWrap.IsAcyclicDigraph(GapObj(D))::Bool

@doc raw"""
    is_aperiodic(D::Digraph) -> Bool

Return `true` if the digraph `D` is aperiodic (the greatest common divisor
of the lengths of all directed cycles is 1).

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_aperiodic(d)
false

julia> d2 = digraph([[2, 3], [3], [1]])
Digraph with 3 vertices, 4 edges

julia> is_aperiodic(d2)
true
```
"""
is_aperiodic(D::Digraph) = DigraphWrap.IsAperiodicDigraph(GapObj(D))::Bool

@doc raw"""
    is_bipartite(D::Digraph) -> Bool

Return `true` if the digraph `D` is bipartite (the vertices can be
partitioned into two sets such that every edge goes from one set to the
other).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1], [1]])
Digraph with 3 vertices, 4 edges

julia> is_bipartite(d)
true

julia> d2 = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_bipartite(d2)
false
```
"""
is_bipartite(D::Digraph) = DigraphWrap.IsBipartiteDigraph(GapObj(D))::Bool

@doc raw"""
    is_complete(D::Digraph) -> Bool

Return `true` if the digraph `D` is complete (every ordered pair of distinct
vertices is joined by an edge in both directions).

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_complete(d)
true

julia> d2 = digraph([[2], [1], [1, 2]])
Digraph with 3 vertices, 4 edges

julia> is_complete(d2)
false
```
"""
is_complete(D::Digraph) = DigraphWrap.IsCompleteDigraph(GapObj(D))::Bool

@doc raw"""
    is_complete_bipartite(D::Digraph) -> Bool

Return `true` if the digraph `D` is a complete bipartite digraph.

# Examples
```jldoctest
julia> d = complete_bipartite_digraph(2, 2)
Digraph with 4 vertices, 8 edges

julia> is_complete_bipartite(d)
true
```
"""
is_complete_bipartite(D::Digraph) = DigraphWrap.IsCompleteBipartiteDigraph(GapObj(D))::Bool

@doc raw"""
    is_complete_multipartite(D::Digraph) -> Bool

Return `true` if the digraph `D` is a complete multipartite digraph.

# Examples
```jldoctest
julia> d = complete_multipartite_digraph([2, 2])
Digraph with 4 vertices, 12 edges

julia> is_complete_multipartite(d)
true
```
"""
is_complete_multipartite(D::Digraph) = DigraphWrap.IsCompleteMultipartiteDigraph(GapObj(D))::Bool

@doc raw"""
    is_tournament(D::Digraph) -> Bool

Return `true` if the digraph `D` is a tournament (every ordered pair of
distinct vertices is joined by exactly one edge).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], [1]])
Digraph with 3 vertices, 4 edges

julia> is_tournament(d)
false
```
"""
is_tournament(D::Digraph) = DigraphWrap.IsTournament(GapObj(D))::Bool

@doc raw"""
    is_symmetric(D::Digraph) -> Bool

Return `true` if the digraph `D` is symmetric (for every edge `(u, v)`,
the reverse edge `(v, u)` is also present).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_symmetric(d)
true
```
"""
is_symmetric(D::Digraph) = DigraphWrap.IsSymmetricDigraph(GapObj(D))::Bool

@doc raw"""
    is_antisymmetric(D::Digraph) -> Bool

Return `true` if the digraph `D` is antisymmetric (there is no pair of
opposite edges between distinct vertices).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_antisymmetric(d)
true

julia> d2 = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_antisymmetric(d2)
false
```
"""
is_antisymmetric(D::Digraph) = DigraphWrap.IsAntisymmetricDigraph(GapObj(D))::Bool

@doc raw"""
    is_transitive(D::Digraph) -> Bool

Return `true` if the digraph `D` is transitive (whenever `(u, v)` and
`(v, w)` are edges, then `(u, w)` is also an edge).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> is_transitive(d)
true

julia> d2 = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_transitive(d2)
false
```
"""
is_transitive(D::Digraph) = DigraphWrap.IsTransitiveDigraph(GapObj(D))::Bool

@doc raw"""
    is_reflexive(D::Digraph) -> Bool

Return `true` if the digraph `D` is reflexive (every vertex has a loop).

# Examples
```jldoctest
julia> d = digraph([[1, 2], [1, 2]])
Digraph with 2 vertices, 4 edges

julia> is_reflexive(d)
true
```
"""
is_reflexive(D::Digraph) = DigraphWrap.IsReflexiveDigraph(GapObj(D))::Bool

@doc raw"""
    is_empty(D::Digraph) -> Bool

Return `true` if the digraph `D` has no edges.

# Examples
```jldoctest
julia> d = null_digraph(3)
Digraph with 3 vertices, 0 edges

julia> is_empty(d)
true
```
"""
is_empty(D::Digraph) = DigraphWrap.IsEmptyDigraph(GapObj(D))::Bool

@doc raw"""
    is_nonempty(D::Digraph) -> Bool

Return `true` if the digraph `D` has at least one edge.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_nonempty(d)
true
```
"""
is_nonempty(D::Digraph) = DigraphWrap.IsNonemptyDigraph(GapObj(D))::Bool

@doc raw"""
    is_null(D::Digraph) -> Bool

Return `true` if the digraph `D` is a null digraph (has no vertices).

# Examples
```jldoctest
julia> d = empty_digraph(0)
Digraph with 0 vertices, 0 edges

julia> is_null(d)
true
```
"""
is_null(D::Digraph) = DigraphWrap.IsNullDigraph(GapObj(D))::Bool

@doc raw"""
    is_eulerian(D::Digraph) -> Bool

Return `true` if the digraph `D` is Eulerian (has an Eulerian trail).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_eulerian(d)
true
```
"""
is_eulerian(D::Digraph) = DigraphWrap.IsEulerianDigraph(GapObj(D))::Bool

@doc raw"""
    is_hamiltonian(D::Digraph) -> Bool

Return `true` if the digraph `D` is Hamiltonian (has a Hamiltonian cycle).

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> is_hamiltonian(d)
true
```
"""
is_hamiltonian(D::Digraph) = DigraphWrap.IsHamiltonianDigraph(GapObj(D))::Bool

@doc raw"""
    is_regular(D::Digraph) -> Bool

Return `true` if the digraph `D` is regular (every vertex has the same
out-degree and in-degree).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_regular(d)
true
```
"""
is_regular(D::Digraph) = DigraphWrap.IsRegularDigraph(GapObj(D))::Bool

@doc raw"""
    is_in_regular(D::Digraph) -> Bool

Return `true` if the digraph `D` is in-regular (every vertex has the same
in-degree).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_in_regular(d)
true
```
"""
is_in_regular(D::Digraph) = DigraphWrap.IsInRegularDigraph(GapObj(D))::Bool

@doc raw"""
    is_out_regular(D::Digraph) -> Bool

Return `true` if the digraph `D` is out-regular (every vertex has the same
out-degree).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_out_regular(d)
true
```
"""
is_out_regular(D::Digraph) = DigraphWrap.IsOutRegularDigraph(GapObj(D))::Bool

@doc raw"""
    is_distance_regular(D::Digraph) -> Bool

Return `true` if the digraph `D` is distance-regular.

# Examples
```jldoctest
julia> d = cycle_digraph(5)
Digraph with 5 vertices, 5 edges

julia> is_distance_regular(d)
false
```
"""
is_distance_regular(D::Digraph) = DigraphWrap.IsDistanceRegularDigraph(GapObj(D))::Bool

@doc raw"""
    is_multi(D::Digraph) -> Bool

Return `true` if the digraph `D` is a multidigraph (has at least one pair
of vertices joined by multiple edges).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 2], [1]])
Digraph with 2 vertices, 4 edges

julia> is_multi(d)
true
```
"""
is_multi(D::Digraph) = DigraphWrap.IsMultiDigraph(GapObj(D))::Bool

@doc raw"""
    is_functional(D::Digraph) -> Bool

Return `true` if the digraph `D` is functional (every vertex has exactly
one out-neighbour).

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> is_functional(d)
true
```
"""
is_functional(D::Digraph) = DigraphWrap.IsFunctionalDigraph(GapObj(D))::Bool

@doc raw"""
    is_permutation(D::Digraph) -> Bool

Return `true` if the digraph `D` is a permutation digraph (functional and
every vertex has exactly one in-neighbour).

# Examples
```jldoctest
julia> d = digraph(3, [1, 2, 3], [2, 3, 1])
Digraph with 3 vertices, 3 edges

julia> is_permutation(d)
true
```
"""
is_permutation(D::Digraph) = DigraphWrap.IsPermutationDigraph(GapObj(D))::Bool

@doc raw"""
    is_equivalence(D::Digraph) -> Bool

Return `true` if the digraph `D` is an equivalence digraph (reflexive,
symmetric, and transitive).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
Digraph with 3 vertices, 9 edges

julia> is_equivalence(d)
true
```
"""
is_equivalence(D::Digraph) = DigraphWrap.IsEquivalenceDigraph(GapObj(D))::Bool

@doc raw"""
    is_preorder(D::Digraph) -> Bool

Return `true` if the digraph `D` is a preorder (reflexive and transitive).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> is_preorder(d)
true
```
"""
is_preorder(D::Digraph) = DigraphWrap.IsPreorderDigraph(GapObj(D))::Bool

@doc raw"""
    is_partial_order(D::Digraph) -> Bool

Return `true` if the digraph `D` is a partial order (antisymmetric,
reflexive, and transitive).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> is_partial_order(d)
true
```
"""
is_partial_order(D::Digraph) = DigraphWrap.IsPartialOrderDigraph(GapObj(D))::Bool

@doc raw"""
    is_meet_semilattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a meet semilattice (a partial order
where every pair of elements has a meet).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2], [3]])
Digraph with 3 vertices, 4 edges

julia> is_meet_semilattice(d)
false
```
"""
is_meet_semilattice(D::Digraph) = DigraphWrap.IsMeetSemilatticeDigraph(GapObj(D))::Bool

@doc raw"""
    is_upper_semimodular(D::Digraph) -> Bool

Return `true` if the digraph `D` is upper semimodular.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 8 edges

julia> is_upper_semimodular(d)
true
```
"""
is_upper_semimodular(D::Digraph) = DigraphWrap.IsUpperSemimodularDigraph(GapObj(D))::Bool

@doc raw"""
    is_distributive_lattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a distributive lattice.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 8 edges

julia> is_distributive_lattice(d)
true
```
"""
is_distributive_lattice(D::Digraph) = DigraphWrap.IsDistributiveLatticeDigraph(GapObj(D))::Bool

@doc raw"""
    is_modular_lattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a modular lattice.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 8 edges

julia> is_modular_lattice(d)
true
```
"""
is_modular_lattice(D::Digraph) = DigraphWrap.IsModularLatticeDigraph(GapObj(D))::Bool

@doc raw"""
    is_biconnected(D::Digraph) -> Bool

Return `true` if the digraph `D` is biconnected (is connected and has
no articulation vertices).

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_biconnected(d)
true
```
"""
is_biconnected(D::Digraph) = DigraphWrap.IsBiconnectedDigraph(GapObj(D))::Bool

@doc raw"""
    is_bridgeless(D::Digraph) -> Bool

Return `true` if the digraph `D` is bridgeless (has no bridges).

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_bridgeless(d)
true
```
"""
is_bridgeless(D::Digraph) = DigraphWrap.IsBridgelessDigraph(GapObj(D))::Bool

@doc raw"""
    is_chain(D::Digraph) -> Bool

Return `true` if the digraph `D` is a chain digraph (a total order).

# Examples
```jldoctest
julia> d = chain_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_chain(d)
true
```
"""
is_chain(D::Digraph) = DigraphWrap.IsChainDigraph(GapObj(D))::Bool

@doc raw"""
    is_cycle(D::Digraph) -> Bool

Return `true` if the digraph `D` is a directed cycle.

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_cycle(d)
true
```
"""
is_cycle(D::Digraph) = DigraphWrap.IsCycleDigraph(GapObj(D))::Bool

@doc raw"""
    is_directed_tree(D::Digraph) -> Bool

Return `true` if the digraph `D` is a directed tree.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_directed_tree(d)
true
```
"""
is_directed_tree(D::Digraph) = DigraphWrap.IsDirectedTree(GapObj(D))::Bool

@doc raw"""
    is_undirected_tree(D::Digraph) -> Bool

Return `true` if the underlying undirected graph of `D` is a tree.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_undirected_tree(d)
false
```
"""
is_undirected_tree(D::Digraph) = DigraphWrap.IsUndirectedTree(GapObj(D))::Bool

@doc raw"""
    is_cayley(D::Digraph) -> Bool

Return `true` if the digraph `D` is a Cayley digraph of some group.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3))
Digraph with 6 vertices, 36 edges

julia> is_cayley(d)
true
```
"""
is_cayley(D::Digraph) = DigraphWrap.IsCayleyDigraph(GapObj(D))::Bool

@doc raw"""
    is_edge_transitive(D::Digraph) -> Bool

Return `true` if the automorphism group of `D` acts transitively on the
edges of `D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_edge_transitive(d)
true
```
"""
is_edge_transitive(D::Digraph) = DigraphWrap.IsEdgeTransitive(GapObj(D))::Bool

@doc raw"""
    is_vertex_transitive(D::Digraph) -> Bool

Return `true` if the automorphism group of `D` acts transitively on the
vertices of `D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_vertex_transitive(d)
true
```
"""
is_vertex_transitive(D::Digraph) = DigraphWrap.IsVertexTransitive(GapObj(D))::Bool

@doc raw"""
    is_2_edge_transitive(D::Digraph) -> Bool

Return `true` if the digraph `D` is 2-edge transitive. Requires the
digraph to have no multiple edges.

# Examples
```jldoctest
julia> d = complete_digraph(4)
Digraph with 4 vertices, 12 edges

julia> is_2_edge_transitive(d)
true
```
"""
is_2_edge_transitive(D::Digraph) = DigraphWrap.Is2EdgeTransitive(GapObj(D))::Bool

@doc raw"""
    is_digraph_core(D::Digraph) -> Bool

Return `true` if the digraph `D` is a core (every endomorphism on `D` is
an automorphism).

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_digraph_core(d)
true
```
"""
is_digraph_core(D::Digraph) = DigraphWrap.IsDigraphCore(GapObj(D))::Bool

@doc raw"""
    is_planar(D::Digraph) -> Bool

Return `true` if the digraph `D` is planar (can be drawn on the plane
without edge crossings).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_planar(d)
true

julia> d2 = complete_digraph(5)
Digraph with 5 vertices, 20 edges

julia> is_planar(d2)
false
```
"""
is_planar(D::Digraph) = DigraphWrap.IsPlanarDigraph(GapObj(D))::Bool

@doc raw"""
    is_outer_planar(D::Digraph) -> Bool

Return `true` if the digraph `D` is outer-planar (can be drawn on the
plane with all vertices on the outer face and no edge crossings).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_outer_planar(d)
true
```
"""
is_outer_planar(D::Digraph) = DigraphWrap.IsOuterPlanarDigraph(GapObj(D))::Bool


@doc raw"""
    is_directed_forest(D::Digraph) -> Bool

Return 	rue if the digraph D is a directed forest (a disjoint union of
directed trees).

# Examples
`jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_directed_forest(d)
true
`
"""
is_directed_forest(D::Digraph) = DigraphWrap.IsDirectedForest(GapObj(D))::Bool

@doc raw"""
    is_undirected_forest(D::Digraph) -> Bool

Return 	rue if the underlying undirected graph of D is a forest.

# Examples
`jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_undirected_forest(d)
true
`
"""
is_undirected_forest(D::Digraph) = DigraphWrap.IsUndirectedForest(GapObj(D))::Bool

@doc raw"""
    is_join_semilattice(D::Digraph) -> Bool

Return 	rue if the digraph D is a join semilattice.

# Examples
`jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_join_semilattice(d)
true
`
"""
is_join_semilattice(D::Digraph) = DigraphWrap.IsJoinSemilatticeDigraph(GapObj(D))::Bool

@doc raw"""
    is_lattice(D::Digraph) -> Bool

Return 	rue if the digraph D is a lattice (both a meet semilattice and
a join semilattice).

# Examples
`jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_lattice(d)
true
`
"""
function is_lattice(D::Digraph)
  return is_meet_semilattice(D) && is_join_semilattice(D)
end

@doc raw"""
    is_cograph(D::Digraph) -> Bool

Return 	rue if the digraph D is a cograph (has no induced subdigraph
isomorphic to a 4-vertex path).

# Examples
`jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_cograph(d)
true
`
"""
is_cograph(D::Digraph) = DigraphWrap.IsCograph(GapObj(D))::Bool

@doc raw"""
    is_lower_semimodular(D::Digraph) -> Bool

Return 	rue if the digraph D represents a lower semimodular lattice.

# Examples
`jldoctest
julia> d = chain_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_lower_semimodular(d)
true
`
"""
is_lower_semimodular(D::Digraph) = DigraphWrap.IsLowerSemimodularDigraph(GapObj(D))::Bool

@doc raw"""
    is_negative_edge_weighted(D::Digraph) -> Bool

Return 	rue if the edge-weighted digraph D has any negative edge weights.

# Examples
`jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [-10]])
Digraph with 2 vertices, 2 edges

julia> is_negative_edge_weighted(d)
true
`
"""
is_negative_edge_weighted(D::Digraph) = DigraphWrap.IsNegativeEdgeWeightedDigraph(GapObj(D))::Bool
