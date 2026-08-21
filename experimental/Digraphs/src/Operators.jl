# ===========================================================================
# Chapter 4: Operators
#
# Wrappers for the operators and operations of Chapter 4 ("Operators") of
# the GAP Digraphs manual:
#
#   4.1 Operators for digraphs
#     - digraph1 = digraph2
#     - digraph1 < digraph2
#     4.1-1 IsSubdigraph(super, sub)
#     4.1-2 IsUndirectedSpanningTree(super, sub)
#     4.1-2 IsUndirectedSpanningForest(super, sub)
#
# The operators == and < delegate to the corresponding GAP operators on the
# underlying GAP objects. These are implemented in the C code of the Digraphs
# package and hence behave exactly as documented: vertex labels are ignored,
# and multiple edges are taken into account.
# ===========================================================================

@doc raw"""
    ==(D1::Digraph, D2::Digraph) -> Bool

Return `true` if the digraphs `D1` and `D2` have the same number of vertices
and the same collection of edges, up to re-ordering of the edge lists.
Multiple edges are taken into account. The vertex labels of `D1` and `D2`
are ignored.

# Examples
```jldoctest
julia> d1 = digraph([[2, 3], [1], [2, 3]])
Digraph with 3 vertices, 5 edges

julia> d2 = digraph([[2, 3], [1], [2, 3]])
Digraph with 3 vertices, 5 edges

julia> d1 == d2
true

julia> d1 == digraph([[2, 3], [], [2]])
false
```
"""
function Base.:(==)(D1::Digraph, D2::Digraph)
  return GapObj(D1) == GapObj(D2)
end

@doc raw"""
    <(D1::Digraph, D2::Digraph) -> Bool

Return `true` if `D1` is less than `D2` in the total order on digraphs used
by GAP. The order is determined as follows:

1. by the number of vertices;
2. if the numbers of vertices are equal, by the number of edges;
3. if the numbers of edges are equal, by comparing the sorted lists of
   edges lexicographically.

# Examples
```jldoctest
julia> complete_digraph(3) < complete_digraph(4)
true

julia> cycle_digraph(4) < complete_digraph(4)
true

julia> digraph([[1], [2]]) < digraph([[2], [1]])
true
```
"""
function Base.:<(D1::Digraph, D2::Digraph)
  return GapObj(D1) < GapObj(D2)
end

@doc raw"""
    is_subdigraph(super::Digraph, sub::Digraph) -> Bool

Return `true` if `sub` is a subdigraph of `super`, and `false` otherwise.

A digraph `sub` is a subdigraph of a digraph `super` if `sub` and `super`
share the same number of vertices, and the collection of edges of `super`
(including repeats) contains the collection of edges of `sub` (including
repeats).

# Examples
```jldoctest
julia> g = digraph([[2, 3], [1], [2, 3]])
Digraph with 3 vertices, 5 edges

julia> h = digraph([[2, 3], [], [2]])
Digraph with 3 vertices, 3 edges

julia> is_subdigraph(g, h)
true

julia> is_subdigraph(h, g)
false

julia> is_subdigraph(complete_digraph(4), cycle_digraph(4))
true

julia> is_subdigraph(cycle_digraph(4), chain_digraph(4))
true

julia> g = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> h = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_subdigraph(g, h)
true

julia> is_subdigraph(h, g)
false
```
"""
function is_subdigraph(super::Digraph, sub::Digraph)
  return DigraphWrap.IsSubdigraph(GapObj(super), GapObj(sub))::Bool
end

@doc raw"""
    is_undirected_spanning_tree(super::Digraph, sub::Digraph) -> Bool

Return `true` if `sub` is an undirected spanning tree of `super`, i.e. `sub`
is a subdigraph of `super` that is an undirected tree. Note that a digraph
whose maximal symmetric subdigraph is not connected has no undirected
spanning trees.

# Examples
```jldoctest
julia> D = complete_digraph(4)
Digraph with 4 vertices, 12 edges

julia> tree = digraph([[3], [4], [1, 4], [2, 3]])
Digraph with 4 vertices, 6 edges

julia> is_undirected_spanning_tree(D, tree)
true
```
"""
function is_undirected_spanning_tree(super::Digraph, sub::Digraph)
  return DigraphWrap.IsUndirectedSpanningTree(GapObj(super), GapObj(sub))::Bool
end

@doc raw"""
    is_undirected_spanning_forest(super::Digraph, sub::Digraph) -> Bool

Return `true` if `sub` is an undirected spanning forest of `super`, i.e.
`sub` is a subdigraph of `super` that is an undirected forest and is not
contained in any larger such subdigraph of `super`. Equivalently, `sub` is
a subdigraph of `super` whose connected components coincide with those of
the maximal symmetric subdigraph of `super`.

# Examples
```jldoctest
julia> D = complete_digraph(4)
Digraph with 4 vertices, 12 edges

julia> forest = empty_digraph(4)
Digraph with 4 vertices, 0 edges

julia> is_undirected_spanning_forest(D, forest)
false

julia> D = digraph_disjoint_union(cycle_digraph(2), cycle_digraph(2))
Digraph with 4 vertices, 4 edges

julia> is_undirected_forest(D) && is_undirected_spanning_forest(D, D)
true
```
"""
function is_undirected_spanning_forest(super::Digraph, sub::Digraph)
  return DigraphWrap.IsUndirectedSpanningForest(GapObj(super), GapObj(sub))::Bool
end

# Make `==`-equal digraphs behave consistently in `Set` and `Dict`: the
# default (identity-based) `isequal` and `hash` for mutable structs would
# otherwise contradict the content-based equality defined above.
Base.isequal(D1::Digraph, D2::Digraph) = D1 == D2

function Base.hash(D::Digraph, h::UInt)
  return hash((nv(D), sort!(digraph_edges(D))), h)
end

