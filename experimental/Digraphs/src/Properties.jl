# ============================================================================
# Properties.jl - full wrapper for the GAP Digraphs package, Chapter 6
# Reference: https://docs.gap-system.org/pkg/digraphs/doc/chap6_mj.html
#
# This file follows the exact order of the GAP manual chapter "Properties of
# digraphs" (6.1 -> 6.2 -> ... -> 6.8) and wraps every version of every
# operation, including GAP aliases (e.g. IsNullDigraph/IsEmptyDigraph,
# IsQuasiorderDigraph/IsPreorderDigraph) and the mutable-copy variants
# (e.g. EdgeWeightsMutableCopy).
# ============================================================================

# ############################################################################
# Section 6.1  Vertex properties
# ############################################################################

# ----------------------------------------------------------------------------
# 6.1-1  DigraphHasAVertex
# ----------------------------------------------------------------------------
@doc raw"""
    has_a_vertex(D::Digraph) -> Bool

Return `true` if the digraph `D` has at least one vertex, and `false`
otherwise.

# Examples
```jldoctest
julia> d = empty_digraph(0)
Digraph with 0 vertices, 0 edges

julia> has_a_vertex(d)
false

julia> d2 = digraph([[]])
Digraph with 1 vertex, 0 edges

julia> has_a_vertex(d2)
true
```
"""
has_a_vertex(D::Digraph) = DigraphWrap.DigraphHasAVertex(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.1-2  DigraphHasNoVertices
# ----------------------------------------------------------------------------
@doc raw"""
    has_no_vertices(D::Digraph) -> Bool

Return `true` if the digraph `D` is the unique digraph with zero vertices,
and `false` otherwise.

# Examples
```jldoctest
julia> d = empty_digraph(0)
Digraph with 0 vertices, 0 edges

julia> has_no_vertices(d)
true
```
"""
has_no_vertices(D::Digraph) = DigraphWrap.DigraphHasNoVertices(GapObj(D))::Bool

# ############################################################################
# Section 6.2  Edge properties
# ############################################################################

# ----------------------------------------------------------------------------
# 6.2-1  DigraphHasLoops
# ----------------------------------------------------------------------------
@doc raw"""
    has_loops(D::Digraph) -> Bool

Return `true` if the digraph `D` has any loops (edges from a vertex to
itself), and `false` otherwise.

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

# ----------------------------------------------------------------------------
# 6.2-2  IsAntiSymmetricDigraph / IsAntisymmetricDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_antisymmetric(D::Digraph) -> Bool

Return `true` if the digraph `D` is antisymmetric, and `false` otherwise.
A digraph is antisymmetric if whenever there is an edge with source `u` and
range `v`, and an edge with source `v` and range `u`, then the vertices `u`
and `v` are equal.

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

# ----------------------------------------------------------------------------
# 6.2-3  IsBipartiteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_bipartite(D::Digraph) -> Bool

Return `true` if the digraph `D` is bipartite, and `false` otherwise. A
digraph is bipartite if and only if the vertices can be partitioned into two
non-empty sets such that the source and range of every edge lie in distinct
sets; equivalently, it is 2-colorable.

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

# ----------------------------------------------------------------------------
# 6.2-4  IsCompleteBipartiteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_complete_bipartite(D::Digraph) -> Bool

Return `true` if the digraph `D` is a complete bipartite digraph, and `false`
otherwise. A digraph is a complete bipartite digraph if it is bipartite and
there is a unique edge with source `i` and range `j` if and only if `i` and
`j` lie in different bicomponents of `D`.

# Examples
```jldoctest
julia> d = complete_bipartite_digraph(2, 2)
Digraph with 4 vertices, 8 edges

julia> is_complete_bipartite(d)
true
```
"""
is_complete_bipartite(D::Digraph) = DigraphWrap.IsCompleteBipartiteDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-5  IsCompleteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_complete(D::Digraph) -> Bool

Return `true` if the digraph `D` is complete, and `false` otherwise. A
digraph is complete if it has no loops and for all distinct vertices `i` and
`j` there is exactly one edge with source `i` and range `j`.

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

# ----------------------------------------------------------------------------
# 6.2-6  IsCompleteMultipartiteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_complete_multipartite(D::Digraph) -> Bool

Return `true` if the digraph `D` is a complete multipartite digraph, and
`false` otherwise. A digraph is a complete multipartite digraph if and only
if its vertices can be partitioned into at least two maximal independent
sets, where every possible edge between these independent sets occurs in the
digraph exactly once.

# Examples
```jldoctest
julia> d = complete_multipartite_digraph([2, 2])
Digraph with 4 vertices, 8 edges

julia> is_complete_multipartite(d)
true
```
"""
is_complete_multipartite(D::Digraph) = DigraphWrap.IsCompleteMultipartiteDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-7  IsEmptyDigraph / IsNullDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_empty(D::Digraph) -> Bool

Return `true` if the digraph `D` is empty (has no edges), and `false`
otherwise.

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
    is_null(D::Digraph) -> Bool

Return `true` if the digraph `D` is a null digraph (has no edges), and
`false` otherwise. `is_null` is a synonym for [`is_empty`](@ref).

# Examples
```jldoctest
julia> d = digraph([[], []])
Digraph with 2 vertices, 0 edges

julia> is_null(d)
true

julia> d2 = digraph([[], [1]])
Digraph with 2 vertices, 1 edge

julia> is_null(d2)
false
```
"""
is_null(D::Digraph) = DigraphWrap.IsNullDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-8  IsEquivalenceDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_equivalence(D::Digraph) -> Bool

Return `true` if the digraph `D` is an equivalence digraph, and `false`
otherwise. A digraph is an equivalence digraph if and only if it satisfies
`is_reflexive`, `is_symmetric` and `is_transitive`.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
Digraph with 3 vertices, 9 edges

julia> is_equivalence(d)
true
```
"""
is_equivalence(D::Digraph) = DigraphWrap.IsEquivalenceDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-9  IsFunctionalDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_functional(D::Digraph) -> Bool

Return `true` if the digraph `D` is functional, and `false` otherwise. A
digraph is functional if every vertex is the source of exactly one edge.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> is_functional(d)
true

julia> d2 = digraph([[1, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> is_functional(d2)
false
```
"""
is_functional(D::Digraph) = DigraphWrap.IsFunctionalDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-10  IsPermutationDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_permutation(D::Digraph) -> Bool

Return `true` if the digraph `D` is a permutation digraph, and `false`
otherwise. A digraph is a permutation digraph if it is functional and every
vertex has exactly one in-neighbour.

# Examples
```jldoctest
julia> d = digraph(3, [1, 2, 3], [2, 3, 1])
Digraph with 3 vertices, 3 edges

julia> is_permutation(d)
true
```
"""
is_permutation(D::Digraph) = DigraphWrap.IsPermutationDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-11  IsMultiDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_multi(D::Digraph) -> Bool

Return `true` if the digraph `D` is a multidigraph, and `false` otherwise. A
multidigraph is one that has at least two edges with equal source and range.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 2], [1]])
Digraph with 2 vertices, 4 edges

julia> is_multi(d)
true
```
"""
is_multi(D::Digraph) = DigraphWrap.IsMultiDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-12  IsNonemptyDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_nonempty(D::Digraph) -> Bool

Return `true` if the digraph `D` is nonempty (has at least one edge), and
`false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_nonempty(d)
true
```
"""
is_nonempty(D::Digraph) = DigraphWrap.IsNonemptyDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-13  IsReflexiveDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_reflexive(D::Digraph) -> Bool

Return `true` if the digraph `D` is reflexive, and `false` otherwise. A
digraph is reflexive if it has a loop at every vertex.

# Examples
```jldoctest
julia> d = digraph([[1, 2], [1, 2]])
Digraph with 2 vertices, 4 edges

julia> is_reflexive(d)
true
```
"""
is_reflexive(D::Digraph) = DigraphWrap.IsReflexiveDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-14  IsSymmetricDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_symmetric(D::Digraph) -> Bool

Return `true` if the digraph `D` is symmetric, and `false` otherwise. A
symmetric digraph is one where for each non-loop edge with source `u` and
range `v` there is a corresponding edge with source `v` and range `u`; if
there are `n` edges with source `u` and range `v`, then there must be
precisely `n` edges with source `v` and range `u`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_symmetric(d)
true
```
"""
is_symmetric(D::Digraph) = DigraphWrap.IsSymmetricDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-15  IsTournament
# ----------------------------------------------------------------------------
@doc raw"""
    is_tournament(D::Digraph) -> Bool

Return `true` if the digraph `D` is a tournament, and `false` otherwise. A
tournament is a digraph which has a unique directed edge (of some
orientation) between any pair of distinct vertices, and no loops.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> is_tournament(d)
true

julia> d2 = digraph([[2, 3], [3], [1]])
Digraph with 3 vertices, 4 edges

julia> is_tournament(d2)
false
```
"""
is_tournament(D::Digraph) = DigraphWrap.IsTournament(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.2-16  IsTransitiveDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_transitive(D::Digraph) -> Bool

Return `true` if the digraph `D` is transitive, and `false` otherwise. A
digraph is transitive if whenever `[i, j]` and `[j, k]` are edges of the
digraph, then `[i, k]` is also an edge.

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

# ############################################################################
# Section 6.3  Edge weights
# ############################################################################

# ----------------------------------------------------------------------------
# 6.3-1  EdgeWeights / EdgeWeightsMutableCopy
# ----------------------------------------------------------------------------
@doc raw"""
    edge_weights(D::Digraph) -> Vector{Vector}

Return the list of lists of edge weights of the edge-weighted digraph `D`,
where `edge_weights(D)[i][j]` is the weight of the `j`-th edge from vertex
`i`, according to the ordering of edges given by `out_neighbours(D)[i]`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [10]])
Digraph with 2 vertices, 2 edges

julia> edge_weights(d)
2-element Vector{Vector}:
 Any[5]
 Any[10]
```
"""
function edge_weights(D::Digraph)
  return Vector{Vector}(DigraphWrap.EdgeWeights(GapObj(D)))
end

@doc raw"""
    edge_weights_mutable_copy(D::Digraph) -> Vector{Vector}

Return a mutable copy of the edge weights of the edge-weighted digraph `D`;
see [`edge_weights`](@ref).

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [10]])
Digraph with 2 vertices, 2 edges

julia> w = edge_weights_mutable_copy(d)
2-element Vector{Vector}:
 Any[5]
 Any[10]
```
"""
function edge_weights_mutable_copy(D::Digraph)
  return Vector{Vector}(DigraphWrap.EdgeWeightsMutableCopy(GapObj(D)))
end

# ----------------------------------------------------------------------------
# 6.3-2  EdgeWeightedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    edge_weighted_digraph(D::Digraph, weights::Vector{<:AbstractVector}) -> Digraph
    edge_weighted_digraph(adj::Vector{<:AbstractVector}, weights::Vector{<:AbstractVector}; mut::Bool=false) -> Digraph

Return an edge-weighted digraph: the digraph `D` (or the digraph with
adjacency list `adj`) together with the edge `weights`, where the `i`-th
entry of `weights` lists the weights of the edges of vertex `i` in the order
of its out-neighbours.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [10]])
Digraph with 2 vertices, 2 edges

julia> edge_weights(d)
2-element Vector{Vector}:
 Any[5]
 Any[10]
```
"""
function edge_weighted_digraph(D::Digraph, weights::Vector{<:AbstractVector})
  return Digraph(DigraphWrap.EdgeWeightedDigraph(GapObj(D),
                                                  GapObj(weights, recursive=true)))
end

function edge_weighted_digraph(adj::Vector{<:AbstractVector},
                               weights::Vector{<:AbstractVector}; mut::Bool=false)
  d = digraph(adj; mut=mut)
  return edge_weighted_digraph(d, weights)
end

# ----------------------------------------------------------------------------
# 6.3-3  EdgeWeightedDigraphTotalWeight
# ----------------------------------------------------------------------------
@doc raw"""
    edge_weighted_digraph_total_weight(D::Digraph) -> GapObj

Return the sum of the weights of the edges of the edge-weighted digraph `D`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [10]])
Digraph with 2 vertices, 2 edges

julia> edge_weighted_digraph_total_weight(d)
15
```
"""
function edge_weighted_digraph_total_weight(D::Digraph)
  return DigraphWrap.EdgeWeightedDigraphTotalWeight(GapObj(D))
end

# ----------------------------------------------------------------------------
# 6.3-4  EdgeWeightedDigraphMinimumSpanningTree
# ----------------------------------------------------------------------------
@doc raw"""
    edge_weighted_digraph_minimum_spanning_tree(D::Digraph) -> Digraph

Return a digraph which is a minimum spanning tree of the edge-weighted
digraph `D`: a subdigraph with the same vertices and a subset of its edges
that form an undirected tree, with the smallest possible total weight among
all spanning trees of `D`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [10]])
Digraph with 2 vertices, 2 edges

julia> edge_weighted_digraph_minimum_spanning_tree(d)
Digraph with 2 vertices, 1 edge
```
"""
function edge_weighted_digraph_minimum_spanning_tree(D::Digraph)
  return Digraph(DigraphWrap.EdgeWeightedDigraphMinimumSpanningTree(GapObj(D)))
end

# ----------------------------------------------------------------------------
# 6.3-5  EdgeWeightedDigraphShortestPaths
# ----------------------------------------------------------------------------
@doc raw"""
    edge_weighted_digraph_shortest_paths(D::Digraph) -> GapObj

Return a GAP record describing the paths of lowest total weight (the
shortest paths) connecting each pair of vertices of the edge-weighted
digraph `D`. The record contains the components `distances`, `parents` and
`edges`, each a list of lists, where the `i`-th entry corresponds to paths
starting at vertex `i`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [3], []], [[5], [10], []])
Digraph with 3 vertices, 2 edges

julia> edge_weighted_digraph_shortest_paths(d)
GAP: rec( distances := [ [ 0, 5, 15 ], [ fail, 0, 10 ], [ fail, fail, 0 ] ], edges := [ [ fail, 1, 1 ], [ fail, fail, 1 ], [ fail, fail, fail ] ], parents := [ [ fail, 1, 1 ], [ fail, fail, 2 ], [ fail, fail, fail ] ] )
```
"""
function edge_weighted_digraph_shortest_paths(D::Digraph)
  return DigraphWrap.EdgeWeightedDigraphShortestPaths(GapObj(D))
end

@doc raw"""
    edge_weighted_digraph_shortest_paths(D::Digraph, source::Int) -> GapObj

Return a GAP record describing the shortest paths from the vertex `source`
to every other vertex of the edge-weighted digraph `D`. The record contains
the components `distances`, `parents` and `edges`, each a list with one
entry for every vertex.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2, 3], [4], [4], []], [[5, 1], [6], [11], []])
Digraph with 4 vertices, 4 edges

julia> edge_weighted_digraph_shortest_paths(d, 1)
GAP: rec( distances := [ 0, 5, 1, 11 ], edges := [ fail, 1, 2, 1 ], parents := [ fail, 1, 1, 2 ] )
```
"""
function edge_weighted_digraph_shortest_paths(D::Digraph, source::Int)
  return DigraphWrap.EdgeWeightedDigraphShortestPaths(GapObj(D), source)
end

# ----------------------------------------------------------------------------
# 6.3-6  EdgeWeightedDigraphShortestPath
# ----------------------------------------------------------------------------
@doc raw"""
    edge_weighted_digraph_shortest_path(D::Digraph, s::Int, t::Int) -> GapObj

Return a directed path from vertex `s` to vertex `t` of the edge-weighted
digraph `D` with the smallest possible total weight, or `fail` if `s == t`
or no such path exists. The output is a pair of lists of the form described
by [`digraph_path`](@ref).

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [3], []], [[5], [10], []])
Digraph with 3 vertices, 2 edges

julia> edge_weighted_digraph_shortest_path(d, 1, 3)
GAP: [ [ 1, 2, 3 ], [ 1, 1 ] ]
```
"""
function edge_weighted_digraph_shortest_path(D::Digraph, s::Int, t::Int)
  return DigraphWrap.EdgeWeightedDigraphShortestPath(GapObj(D), s, t)
end

# ----------------------------------------------------------------------------
# 6.3-7  DigraphMaximumFlow
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_maximum_flow(D::Digraph, s::Int, t::Int) -> GapObj

Return the maximum flow from vertex `s` to vertex `t` in the edge-weighted
digraph `D` as a list of lists, where each entry is the flow on the edge in
the corresponding position of `out_neighbours(D)`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2, 2], [3], []], [[3, 2], [1], []])
Digraph with 3 vertices, 3 edges

julia> digraph_maximum_flow(d, 1, 3)
GAP: [ [ 1, 0 ], [ 1 ], [  ] ]
```
"""
function digraph_maximum_flow(D::Digraph, s::Int, t::Int)
  return DigraphWrap.DigraphMaximumFlow(GapObj(D), s, t)
end

# ----------------------------------------------------------------------------
# 6.3-8  DigraphMinimumCut
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_minimum_cut(D::Digraph, s::Int, t::Int) -> Vector{Vector{Int}}

Return a list of two lists representing the components of the minimum
`s`-`t` cut of the edge-weighted digraph `D` with distinct vertices `s` and
`t`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2, 2], [3], []], [[3, 2], [1], []])
Digraph with 3 vertices, 3 edges

julia> digraph_minimum_cut(d, 1, 3)
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3]
```
"""
function digraph_minimum_cut(D::Digraph, s::Int, t::Int)
  return Vector{Vector{Int}}(DigraphWrap.DigraphMinimumCut(GapObj(D), s, t))
end

# ----------------------------------------------------------------------------
# 6.3-9  DigraphMinimumCutSet
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_minimum_cut_set(D::Digraph, s::Int, t::Int) -> Vector{Tuple{Int,Int}}

Return the set of edges (as a list of source-range pairs) of the minimum
`s`-`t` cut of the edge-weighted digraph `D` with distinct vertices `s` and
`t`.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2, 2], [3], []], [[3, 2], [1], []])
Digraph with 3 vertices, 3 edges

julia> digraph_minimum_cut_set(d, 1, 3)
1-element Vector{Tuple{Int64, Int64}}:
 (2, 3)
```
"""
function digraph_minimum_cut_set(D::Digraph, s::Int, t::Int)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphMinimumCutSet(GapObj(D), s, t))
end

# ----------------------------------------------------------------------------
# 6.3-10  RandomUniqueEdgeWeightedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    random_unique_edge_weighted_digraph(n::Integer; mut::Bool=false) -> Digraph
    random_unique_edge_weighted_digraph(n::Integer, p::AbstractFloat; mut::Bool=false) -> Digraph

Return a random edge-weighted digraph with `n` vertices, whose edge weights
are the random unique weights from the set `1:m`, where `m` is the number of
edges of the digraph. If `p` is given, an edge exists between each pair of
vertices with probability approximately `p`; otherwise a random probability
is assumed. The keyword `mut` selects the GAP filter for the created
digraph. See [`random_digraph`](@ref) for more details.

# Examples
```jldoctest
julia> d = random_unique_edge_weighted_digraph(4, 0.0)
Digraph with 4 vertices, 0 edges

julia> edge_weighted_digraph_total_weight(d)
0
```
"""
function random_unique_edge_weighted_digraph(n::Integer; mut::Bool=false)
  return Digraph(DigraphWrap.RandomUniqueEdgeWeightedDigraph(_filt(mut), Int(n)))
end

function random_unique_edge_weighted_digraph(n::Integer, p::AbstractFloat; mut::Bool=false)
  return Digraph(DigraphWrap.RandomUniqueEdgeWeightedDigraph(_filt(mut), Int(n), GapObj(p)))
end

# ############################################################################
# Section 6.4  Orders
# ############################################################################

# ----------------------------------------------------------------------------
# 6.4-1  IsPreorderDigraph / IsQuasiorderDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_preorder(D::Digraph) -> Bool

Return `true` if the digraph `D` is a preorder digraph (reflexive and
transitive), and `false` otherwise. A preorder digraph corresponds to the
preorder relation `x <= y` if and only if `[x, y]` is an edge of `D`.

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
    is_quasiorder(D::Digraph) -> Bool

Return `true` if the digraph `D` is a quasiorder digraph, and `false`
otherwise. `is_quasiorder` is a synonym for [`is_preorder`](@ref).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> is_quasiorder(d)
true
```
"""
is_quasiorder(D::Digraph) = DigraphWrap.IsQuasiorderDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.4-2  IsPartialOrderDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_partial_order(D::Digraph) -> Bool

Return `true` if the digraph `D` is a partial order digraph (reflexive,
antisymmetric and transitive), and `false` otherwise. A partial order
digraph corresponds to the partial order relation `x <= y` if and only if
`[x, y]` is an edge of `D`.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> is_partial_order(d)
true
```
"""
is_partial_order(D::Digraph) = DigraphWrap.IsPartialOrderDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.4-3  IsMeetSemilatticeDigraph / IsJoinSemilatticeDigraph / IsLatticeDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_meet_semilattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a meet semilattice, and `false`
otherwise. A digraph is a meet semilattice if it is a partial order and
every pair of vertices has a greatest lower bound (meet).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2], [3]])
Digraph with 3 vertices, 5 edges

julia> is_meet_semilattice(d)
true

julia> d2 = digraph([[1, 2], [2], [3]])
Digraph with 3 vertices, 4 edges

julia> is_meet_semilattice(d2)
false
```
"""
is_meet_semilattice(D::Digraph) = DigraphWrap.IsMeetSemilatticeDigraph(GapObj(D))::Bool

@doc raw"""
    is_join_semilattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a join semilattice, and `false`
otherwise. A digraph is a join semilattice if it is a partial order and
every pair of vertices has a least upper bound (join).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_join_semilattice(d)
true

julia> d2 = digraph([[1, 2, 3], [2], [3]])
Digraph with 3 vertices, 5 edges

julia> is_join_semilattice(d2)
false
```
"""
is_join_semilattice(D::Digraph) = DigraphWrap.IsJoinSemilatticeDigraph(GapObj(D))::Bool

@doc raw"""
    is_lattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a lattice (both a meet and a join
semilattice), and `false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_lattice(d)
true
```
"""
function is_lattice(D::Digraph)
  return is_meet_semilattice(D) && is_join_semilattice(D)
end

# ----------------------------------------------------------------------------
# 6.4-4  DigraphMeetTable / DigraphJoinTable
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_meet_table(D::Digraph) -> Matrix{Int}

Return the meet table of the digraph `D`. Throw an error if `D` is not a
meet semilattice.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> digraph_meet_table(d)
4×4 Matrix{Int64}:
 1  1  1  1
 1  2  1  2
 1  1  3  3
 1  2  3  4
```
"""
function digraph_meet_table(D::Digraph)
  result = DigraphWrap.DigraphMeetTable(GapObj(D))
  @req result != GAP.Globals.fail "digraph_meet_table: the digraph is not a meet semilattice"
  return Matrix{Int}(result)
end

@doc raw"""
    digraph_join_table(D::Digraph) -> Matrix{Int}

Return the join table of the digraph `D`. Throw an error if `D` is not a
join semilattice.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> digraph_join_table(d)
4×4 Matrix{Int64}:
 1  2  3  4
 2  2  4  4
 3  4  3  4
 4  4  4  4
```
"""
function digraph_join_table(D::Digraph)
  result = DigraphWrap.DigraphJoinTable(GapObj(D))
  @req result != GAP.Globals.fail "digraph_join_table: the digraph is not a join semilattice"
  return Matrix{Int}(result)
end

# ----------------------------------------------------------------------------
# 6.4-5  IsOrderIdeal
# ----------------------------------------------------------------------------
@doc raw"""
    is_order_ideal(D::Digraph, subset::Vector{<:Integer}) -> Bool

Return `true` if `subset` is closed under reachability in the partial order
digraph `D`: every vertex of `D` that can be reached from a vertex of
`subset` by following edges of `D` is again contained in `subset`. Throw an
error if `D` is not a partial order digraph.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> is_order_ideal(d, [3])
true

julia> is_order_ideal(d, [1])
false
```
"""
function is_order_ideal(D::Digraph, subset::Vector{<:Integer})
  return DigraphWrap.IsOrderIdeal(GapObj(D), GapObj(collect(Int, subset)))::Bool
end

# ----------------------------------------------------------------------------
# 6.4-6  IsOrderFilter
# ----------------------------------------------------------------------------
@doc raw"""
    is_order_filter(D::Digraph, subset::Vector{<:Integer}) -> Bool

Return `true` if `subset` is an order filter of the partial order digraph
`D`, i.e. if `subset` is an order ideal of the reverse of `D`. Throw an
error if `D` is not a partial order digraph.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> is_order_filter(d, [1])
true

julia> is_order_filter(d, [2])
false
```
"""
function is_order_filter(D::Digraph, subset::Vector{<:Integer})
  return DigraphWrap.IsOrderFilter(GapObj(D), GapObj(collect(Int, subset)))::Bool
end

# ----------------------------------------------------------------------------
# 6.4-7  IsUpperSemimodularDigraph / IsLowerSemimodularDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_upper_semimodular(D::Digraph) -> Bool

Return `true` if the digraph `D` represents an upper semimodular lattice,
and `false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_upper_semimodular(d)
true
```
"""
is_upper_semimodular(D::Digraph) = DigraphWrap.IsUpperSemimodularDigraph(GapObj(D))::Bool

@doc raw"""
    is_lower_semimodular(D::Digraph) -> Bool

Return `true` if the digraph `D` represents a lower semimodular lattice, and
`false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_lower_semimodular(d)
true
```
"""
is_lower_semimodular(D::Digraph) = DigraphWrap.IsLowerSemimodularDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.4-8  IsDistributiveLatticeDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_distributive_lattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a distributive lattice digraph, and
`false` otherwise. A distributive lattice digraph is a lattice digraph in
which the meet and join operations distribute over each other.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_distributive_lattice(d)
true
```
"""
is_distributive_lattice(D::Digraph) = DigraphWrap.IsDistributiveLatticeDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.4-9  IsModularLatticeDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_modular_lattice(D::Digraph) -> Bool

Return `true` if the digraph `D` is a modular lattice digraph, and `false`
otherwise. A modular lattice digraph is a lattice digraph in which the
5-element lattice `N5` is not embeddable as a lattice.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> is_modular_lattice(d)
true
```
"""
is_modular_lattice(D::Digraph) = DigraphWrap.IsModularLatticeDigraph(GapObj(D))::Bool

# ############################################################################
# Section 6.5  Regularity
# ############################################################################

# ----------------------------------------------------------------------------
# 6.5-1  IsInRegularDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_in_regular(D::Digraph) -> Bool

Return `true` if there is an integer `n` such that every vertex of `D` has
exactly `n` edges terminating in it, and `false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_in_regular(d)
true
```
"""
is_in_regular(D::Digraph) = DigraphWrap.IsInRegularDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.5-2  IsOutRegularDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_out_regular(D::Digraph) -> Bool

Return `true` if there is an integer `n` such that every vertex of `D` has
exactly `n` edges starting at it, and `false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_out_regular(d)
true
```
"""
is_out_regular(D::Digraph) = DigraphWrap.IsOutRegularDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.5-3  IsRegularDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_regular(D::Digraph) -> Bool

Return `true` if there is an integer `n` such that every vertex of `D` has
exactly `n` edges starting and terminating at it, and `false` otherwise. In
other words, the property is true if `D` is both in-regular and out-regular.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_regular(d)
true
```
"""
is_regular(D::Digraph) = DigraphWrap.IsRegularDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.5-4  IsDistanceRegularDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_distance_regular(D::Digraph) -> Bool

Return `true` if the digraph `D` is distance-regular, and `false` otherwise.
If `D` is not symmetric or not connected, the property is `false`.

# Examples
```jldoctest
julia> d = cycle_digraph(5)
Digraph with 5 vertices, 5 edges

julia> is_distance_regular(d)
false
```
"""
is_distance_regular(D::Digraph) = DigraphWrap.IsDistanceRegularDigraph(GapObj(D))::Bool

# ############################################################################
# Section 6.6  Connectivity and cycles
# ############################################################################

# ----------------------------------------------------------------------------
# 6.6-1  IsAcyclicDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_acyclic(D::Digraph) -> Bool

Return `true` if the digraph `D` is acyclic (every directed cycle on the
digraph is trivial), and `false` otherwise.

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

# ----------------------------------------------------------------------------
# 6.6-2  IsChainDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_chain(D::Digraph) -> Bool

Return `true` if the digraph `D` is a chain digraph (isomorphic to the chain
digraph with the same number of vertices), and `false` otherwise.

# Examples
```jldoctest
julia> d = chain_digraph(3)
Digraph with 3 vertices, 2 edges

julia> is_chain(d)
true
```
"""
is_chain(D::Digraph) = DigraphWrap.IsChainDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-3  IsConnectedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_connected(D::Digraph) -> Bool

Return `true` if the digraph `D` is weakly connected, and `false` otherwise.
A digraph is weakly connected if it is possible to travel from any vertex to
any other vertex by traversing edges in either direction.

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

# ----------------------------------------------------------------------------
# 6.6-4  IsBiconnectedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_biconnected(D::Digraph) -> Bool

Return `true` if the digraph `D` is biconnected, and `false` otherwise. A
connected digraph is biconnected if it is still connected when any vertex is
removed. Multiple edges are ignored by this method.

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_biconnected(d)
true
```
"""
is_biconnected(D::Digraph) = DigraphWrap.IsBiconnectedDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-5  IsBridgelessDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_bridgeless(D::Digraph) -> Bool

Return `true` if the digraph `D` is bridgeless, and `false` otherwise. A
connected digraph is bridgeless if it is still connected when any edge is
removed. Multiple edges are ignored by this method.

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_bridgeless(d)
true
```
"""
is_bridgeless(D::Digraph) = DigraphWrap.IsBridgelessDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-6  IsStronglyConnectedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_strongly_connected(D::Digraph) -> Bool

Return `true` if the digraph `D` is strongly connected, and `false`
otherwise. A digraph is strongly connected if there is a directed path from
every vertex to every other vertex.

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

# ----------------------------------------------------------------------------
# 6.6-7  IsAperiodicDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_aperiodic(D::Digraph) -> Bool

Return `true` if the digraph `D` is aperiodic, i.e. if its period is equal
to 1, and `false` otherwise.

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

# ----------------------------------------------------------------------------
# 6.6-8  IsDirectedTree / IsDirectedForest
# ----------------------------------------------------------------------------
@doc raw"""
    is_directed_tree(D::Digraph) -> Bool

Return `true` if the digraph `D` is a directed tree, and `false` otherwise.
A directed tree is an acyclic digraph with precisely one source, without
multiple edges, and such that no two vertices share an out-neighbour.

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
    is_directed_forest(D::Digraph) -> Bool

Return `true` if the digraph `D` is a directed forest, and `false`
otherwise. A directed forest is a digraph with at least one vertex, each of
whose connected components is a directed tree.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_directed_forest(d)
true

julia> d2 = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_directed_forest(d2)
false
```
"""
is_directed_forest(D::Digraph) = DigraphWrap.IsDirectedForest(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-9  IsUndirectedTree / IsUndirectedForest
# ----------------------------------------------------------------------------
@doc raw"""
    is_undirected_tree(D::Digraph) -> Bool

Return `true` if the digraph `D` is an undirected tree, and `false`
otherwise. An undirected tree is a symmetric digraph without loops, in which
for any pair of distinct vertices `u` and `v`, there is exactly one directed
path from `u` to `v`.

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
    is_undirected_forest(D::Digraph) -> Bool

Return `true` if the digraph `D` is an undirected forest, and `false`
otherwise. An undirected forest is a digraph, each of whose connected
components is an undirected tree.

# Examples
```jldoctest
julia> d = digraph([[], [], []])
Digraph with 3 vertices, 0 edges

julia> is_undirected_forest(d)
true
```
"""
is_undirected_forest(D::Digraph) = DigraphWrap.IsUndirectedForest(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-10  IsEulerianDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_eulerian(D::Digraph) -> Bool

Return `true` if the digraph `D` is Eulerian, and `false` otherwise. A
connected digraph is called Eulerian if there exists a directed circuit on
the digraph which includes every edge exactly once. The empty digraph with
at most one vertex is considered to be Eulerian.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> is_eulerian(d)
true
```
"""
is_eulerian(D::Digraph) = DigraphWrap.IsEulerianDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-11  IsHamiltonianDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_hamiltonian(D::Digraph) -> Bool

Return `true` if the digraph `D` is Hamiltonian, and `false` otherwise. A
digraph with `n` vertices is Hamiltonian if it has a directed cycle of
length `n`. The empty digraphs on 0 and 1 vertices are considered to be
Hamiltonian.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> is_hamiltonian(d)
true
```
"""
is_hamiltonian(D::Digraph) = DigraphWrap.IsHamiltonianDigraph(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.6-12  IsCycleDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_cycle(D::Digraph) -> Bool

Return `true` if the digraph `D` is a cycle digraph (isomorphic to the cycle
digraph with the same number of vertices), and `false` otherwise. A digraph
is a cycle if and only if it is strongly connected and has the same number
of edges as vertices.

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

julia> is_cycle(d)
true
```
"""
is_cycle(D::Digraph) = DigraphWrap.IsCycleDigraph(GapObj(D))::Bool

# ############################################################################
# Section 6.7  Planarity
# ############################################################################

# ----------------------------------------------------------------------------
# 6.7-1  IsPlanarDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_planar(D::Digraph) -> Bool

Return `true` if the digraph `D` is planar, and `false` otherwise. A planar
digraph is a digraph that can be embedded in the plane in such a way that
its edges do not intersect. The directions and multiplicities of the edges
are ignored.

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

# ----------------------------------------------------------------------------
# 6.7-2  IsOuterPlanarDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    is_outer_planar(D::Digraph) -> Bool

Return `true` if the digraph `D` is outer planar, and `false` otherwise. An
outer planar digraph is a digraph that can be embedded in the plane in such
an way that its edges do not intersect and all vertices belong to the
unbounded face of the embedding. The directions and multiplicities of the
edges are ignored.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_outer_planar(d)
true
```
"""
is_outer_planar(D::Digraph) = DigraphWrap.IsOuterPlanarDigraph(GapObj(D))::Bool

# ############################################################################
# Section 6.8  Homomorphisms and transformations
# ############################################################################

# ----------------------------------------------------------------------------
# 6.8-1  IsDigraphCore
# ----------------------------------------------------------------------------
@doc raw"""
    is_digraph_core(D::Digraph) -> Bool

Return `true` if the digraph `D` is a core, and `false` otherwise. A digraph
`D` is a core if and only if every endomorphism on `D` is an automorphism on
`D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_digraph_core(d)
true
```
"""
is_digraph_core(D::Digraph) = DigraphWrap.IsDigraphCore(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.8-2  IsEdgeTransitive
# ----------------------------------------------------------------------------
@doc raw"""
    is_edge_transitive(D::Digraph) -> Bool

Return `true` if the digraph `D` (without multiple edges) is edge
transitive, and `false` otherwise. A digraph is edge transitive if its
automorphism group acts transitively on its edges.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_edge_transitive(d)
true
```
"""
is_edge_transitive(D::Digraph) = DigraphWrap.IsEdgeTransitive(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.8-3  IsVertexTransitive
# ----------------------------------------------------------------------------
@doc raw"""
    is_vertex_transitive(D::Digraph) -> Bool

Return `true` if the digraph `D` is vertex transitive, and `false`
otherwise. A digraph is vertex transitive if its automorphism group acts
transitively on its vertices.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> is_vertex_transitive(d)
true
```
"""
is_vertex_transitive(D::Digraph) = DigraphWrap.IsVertexTransitive(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# 6.8-4  Is2EdgeTransitive
# ----------------------------------------------------------------------------
@doc raw"""
    is_2_edge_transitive(D::Digraph) -> Bool

Return `true` if the digraph `D` (without multiple edges) is 2-edge
transitive, and `false` otherwise. A digraph is 2-edge transitive if its
automorphism group acts transitively on 2-edges: triples `(u, v, w)` of
distinct vertices such that `(u, v)` and `(v, w)` are edges. Throw an error
if `D` has multiple edges.

# Examples
```jldoctest
julia> d = complete_digraph(4)
Digraph with 4 vertices, 12 edges

julia> is_2_edge_transitive(d)
true
```
"""
is_2_edge_transitive(D::Digraph) = DigraphWrap.Is2EdgeTransitive(GapObj(D))::Bool

# ----------------------------------------------------------------------------
# Functions kept from the original wrapper, not part of GAP manual Chapter 6
# ----------------------------------------------------------------------------

@doc raw"""
    is_cayley(D::Digraph) -> Bool

Return `true` if the digraph `D` is a Cayley digraph of some group, and
`false` otherwise.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3))
Digraph with 6 vertices, 12 edges

julia> is_cayley(d)
true
```
"""
is_cayley(D::Digraph) = DigraphWrap.IsCayleyDigraph(GapObj(D))::Bool

@doc raw"""
    is_negative_edge_weighted(D::Digraph) -> Bool

Return `true` if the edge-weighted digraph `D` has any negative edge
weights, and `false` otherwise.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [-10]])
Digraph with 2 vertices, 2 edges

julia> is_negative_edge_weighted(d)
true
```
"""
is_negative_edge_weighted(D::Digraph) = DigraphWrap.IsNegativeEdgeWeightedDigraph(GapObj(D))::Bool
