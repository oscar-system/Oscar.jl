# ===========================================================================
# Chapter 5: Attributes and operations
#
# Wrappers for Chapter 5 ("Attributes and operations") of the GAP Digraphs
# manual, in the order of the manual. Functions that are not part of
# Chapter 5 (e.g. for edge-weighted digraphs) are collected at the end.
# ===========================================================================
@doc raw"""
    nv(D::Digraph) -> Int

Return the number of vertices of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> nv(d)
2
```
"""
function nv(D::Digraph)
  return Int(D.nv)
end

@doc raw"""
    ne(D::Digraph) -> Int

Return the number of edges of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> ne(d)
2
```
"""
function ne(D::Digraph)
  return Int(D.ne)
end

@doc raw"""
    digraph_vertices(D::Digraph) -> Vector{Int}

Return the list of vertices of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> digraph_vertices(d)
2-element Vector{Int64}:
 1
 2
```
"""
function digraph_vertices(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphVertices(GapObj(D)))
end

@doc raw"""
    digraph_edges(D::Digraph) -> Vector{Tuple{Int,Int}}

Return the edges of the digraph `D` as a list of source-target pairs.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> digraph_edges(d)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (2, 1)
```
"""
function digraph_edges(D::Digraph)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphEdges(GapObj(D)))
end

@doc raw"""
    digraph_nr_adjacencies(D::Digraph) -> Int

Return the total number of adjacencies (out-neighbour entries) of `D`,
counting multiplicities.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> digraph_nr_adjacencies(d)
2
```
"""
function digraph_nr_adjacencies(D::Digraph)
  return Int(DigraphWrap.DigraphNrAdjacencies(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_nr_adjacencies_without_loops(D::Digraph) -> Int

Return the number of unordered pairs of distinct vertices connected by
at least one edge in D.

# Examples
```jldoctest
julia> d = digraph([[1, 2], [2]])
Digraph with 2 vertices, 3 edges

julia> digraph_nr_adjacencies_without_loops(d)
1
```
"""
function digraph_nr_adjacencies_without_loops(D::Digraph)
  return Int(DigraphWrap.DigraphNrAdjacenciesWithoutLoops(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_nr_loops(D::Digraph) -> Int

Return the number of loops of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[1, 2], [2]])
Digraph with 2 vertices, 3 edges

julia> digraph_nr_loops(d)
2
```
"""
function digraph_nr_loops(D::Digraph)
  return Int(DigraphWrap.DigraphNrLoops(GapObj(D))::GapInt)
end
@doc raw"""
    digraph_sinks(D::Digraph) -> Vector{Int}

Return the list of sink vertices (vertices with out-degree 0) of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], []])
Digraph with 2 vertices, 1 edge

julia> digraph_sinks(d)
1-element Vector{Int64}:
 2
```
"""
function digraph_sinks(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphSinks(GapObj(D)))
end

@doc raw"""
    digraph_sources(D::Digraph) -> Vector{Int}

Return the list of source vertices (vertices with in-degree 0) of `D`.

# Examples
```jldoctest
julia> d = digraph([[], [1]])
Digraph with 2 vertices, 1 edge

julia> digraph_sources(d)
1-element Vector{Int64}:
 2
```
"""
function digraph_sources(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphSources(GapObj(D)))
end

@doc raw"""
    digraph_topological_sort(D::Digraph) -> Vector{Int}

Return a topological ordering of the vertices of `D`. Throw an error if
`D` contains a directed cycle of length greater than 1.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_topological_sort(d)
3-element Vector{Int64}:
 3
 2
 1
```
"""
function digraph_topological_sort(D::Digraph)
  result = DigraphWrap.DigraphTopologicalSort(GapObj(D))
  @req result != GAP.Globals.fail "digraph_topological_sort: the digraph has no topological sort because it contains a directed cycle of length greater than 1"
  return Vector{Int}(result)
end

function _gap_labels(x)
  if GAP.Globals.IsString(x)
    return String(x)
  elseif GAP.Globals.IsInt(x)
    return Int(x)
  elseif GAP.Globals.IsBool(x)
    return Bool(x)
  elseif GAP.Globals.IsList(x)
    return Any[_gap_labels(x[i]) for i in 1:length(x)]
  else
    return x
  end
end

@doc raw"""
    digraph_vertex_label(D::Digraph, i::Int) -> Any

Return the label of vertex `i` of the digraph `D`. If no label has been set
for vertex `i`, then the default label is `i` itself.

# Examples
```jldoctest
julia> D = digraph(["a", "b", "c"], [], [])
Digraph with 3 vertices, 0 edges

julia> digraph_vertex_label(D, 2)
"b"
```
"""
function digraph_vertex_label(D::Digraph, i::Int)
  return _gap_labels(DigraphWrap.DigraphVertexLabel(GapObj(D), i))
end

@doc raw"""
    set_digraph_vertex_label!(D::Digraph, i::Int, obj) -> Nothing

Set the label of vertex `i` of the digraph `D` to the object `obj`. The label
of a vertex can be changed an arbitrary number of times.

# Examples
```jldoctest
julia> D = digraph(["a", "b", "c"], [], [])
Digraph with 3 vertices, 0 edges

julia> set_digraph_vertex_label!(D, 2, "d")

julia> digraph_vertex_label(D, 2)
"d"
```
"""
function set_digraph_vertex_label!(D::Digraph, i::Int, obj)
  return DigraphWrap.SetDigraphVertexLabel(GapObj(D), i, GapObj(obj, recursive = true))
end

@doc raw"""
    digraph_vertex_labels(D::Digraph) -> Vector{Any}

Return a copy of the labels of the vertices of the digraph `D`. If no labels
have been set, then the default labels are the vertices themselves.

# Examples
```jldoctest
julia> D = digraph(["a", "b", "c"], [], [])
Digraph with 3 vertices, 0 edges

julia> digraph_vertex_labels(D)
3-element Vector{Any}:
 "a"
 "b"
 "c"
```
"""
function digraph_vertex_labels(D::Digraph)
  return Vector{Any}(_gap_labels(DigraphWrap.DigraphVertexLabels(GapObj(D))))
end

@doc raw"""
    set_digraph_vertex_labels!(D::Digraph, list) -> Nothing

Set the labels of the vertices of the digraph `D` to the list `list` of
arbitrary objects, which must have the same length as the number of vertices
of `D`.

# Examples
```jldoctest
julia> D = digraph(["a", "b", "c"], [], [])
Digraph with 3 vertices, 0 edges

julia> set_digraph_vertex_labels!(D, ["h", "i", "j"])

julia> digraph_vertex_labels(D)
3-element Vector{Any}:
 "h"
 "i"
 "j"
```
"""
function set_digraph_vertex_labels!(D::Digraph, list)
  return DigraphWrap.SetDigraphVertexLabels(GapObj(D), GapObj(list, recursive = true))
end

@doc raw"""
    digraph_edge_label(D::Digraph, i::Int, j::Int) -> Any

Return the label of the edge from vertex `i` to vertex `j` of the digraph `D`,
which must have no multiple edges. If no label has been set for the edge,
then the default label is `1`.

# Examples
```jldoctest
julia> D = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> digraph_edge_label(D, 1, 2)
1
```
"""
function digraph_edge_label(D::Digraph, i::Int, j::Int)
  return _gap_labels(DigraphWrap.DigraphEdgeLabel(GapObj(D), i, j))
end

@doc raw"""
    set_digraph_edge_label!(D::Digraph, i::Int, j::Int, obj) -> Nothing

Set the label of the edge from vertex `i` to vertex `j` of the digraph `D`
(which must have no multiple edges) to the object `obj`.

# Examples
```jldoctest
julia> D = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> set_digraph_edge_label!(D, 1, 2, [42])

julia> digraph_edge_label(D, 1, 2)
1-element Vector{Any}:
 42
```
"""
function set_digraph_edge_label!(D::Digraph, i::Int, j::Int, obj)
  return DigraphWrap.SetDigraphEdgeLabel(GapObj(D), i, j, GapObj(obj, recursive = true))
end

@doc raw"""
    digraph_edge_labels(D::Digraph) -> Vector{Vector{Any}}

Return a copy of the labels of the edges of the digraph `D`, which must have
no multiple edges, as a list of lists such that position `[i][j]` holds the
label on the edge from vertex `i` to vertex `out_neighbours(D)[i][j]`. If no
label has been set for an edge, then the default label is `1`.

# Examples
```jldoctest
julia> D = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> digraph_edge_labels(D)
2-element Vector{Vector{Any}}:
 [1]
 [1]
```
"""
function digraph_edge_labels(D::Digraph)
  result = DigraphWrap.DigraphEdgeLabels(GapObj(D))
  return Vector{Vector{Any}}(_gap_labels(result))
end

@doc raw"""
    set_digraph_edge_labels!(D::Digraph, labels) -> Nothing

    set_digraph_edge_labels!(D::Digraph, func::Function) -> Nothing

Set the labels of the edges of the digraph `D` (which must have no multiple
edges) to the list `labels` of lists of arbitrary objects, or to the labels
computed by the binary function `func`, which is passed two vertices `i` and
`j` and returns the label for the edge between them.

# Examples
```jldoctest
julia> D = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> set_digraph_edge_labels!(D, (i, j) -> 7)

julia> digraph_edge_labels(D)
2-element Vector{Vector{Any}}:
 [7]
 [7]
```
"""
function set_digraph_edge_labels!(D::Digraph, labels)
  return DigraphWrap.SetDigraphEdgeLabels(GapObj(D), GapObj(labels, recursive = true))
end

function set_digraph_edge_labels!(D::Digraph, func::Function)
  return DigraphWrap.SetDigraphEdgeLabels(GapObj(D), _gap_wrap_function(func))
end

@doc raw"""
    digraph_in_edges(D::Digraph, v::Int) -> Vector{Tuple{Int, Int}}

Return the list of all edges of the digraph `D` which have vertex `v` as
their range.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> digraph_in_edges(d, 2)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (1, 2)
```
"""
function digraph_in_edges(D::Digraph, v::Int)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphInEdges(GapObj(D), v))
end

@doc raw"""
    digraph_out_edges(D::Digraph, v::Int) -> Vector{Tuple{Int, Int}}

Return the list of all edges of the digraph `D` which have vertex `v` as
their source.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> digraph_out_edges(d, 1)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (1, 2)
```
"""
function digraph_out_edges(D::Digraph, v::Int)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphOutEdges(GapObj(D), v))
end

@doc raw"""
    is_digraph_edge(D::Digraph, e::Vector{<:Integer}) -> Bool

    is_digraph_edge(D::Digraph, u::Int, v::Int) -> Bool

Return `true` if `[u, v]` (or the pair given by the two-element list `e`) is
an edge of the digraph `D`, and `false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [6], [], [3], [], [1]])
Digraph with 6 vertices, 5 edges

julia> is_digraph_edge(d, [1, 2])
true

julia> is_digraph_edge(d, 1, 1)
false
```
"""
function is_digraph_edge(D::Digraph, e::Vector{<:Integer})
  return DigraphWrap.IsDigraphEdge(GapObj(D), e)::Bool
end

function is_digraph_edge(D::Digraph, u::Int, v::Int)
  return DigraphWrap.IsDigraphEdge(GapObj(D), u, v)::Bool
end

@doc raw"""
    is_matching(D::Digraph, edges::Vector) -> Bool

Return `true` if `edges` is a matching of the digraph `D`, i.e. a list of
pairs of vertices no two of which are incident to the same vertex.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_matching(d, [[1, 2]])
true

julia> is_matching(d, [[1, 2], [2, 1]])
false
```
"""
function is_matching(D::Digraph, edges::Vector)
  return DigraphWrap.IsMatching(GapObj(D), GapObj(edges, recursive = true))::Bool
end

@doc raw"""
    is_maximal_matching(D::Digraph, edges::Vector) -> Bool

Return `true` if `edges` is a maximal matching of the digraph `D`, and
`false` otherwise. See [`is_matching`](@ref) for the definition of a
matching.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_maximal_matching(d, [[1, 2]])
true
```
"""
function is_maximal_matching(D::Digraph, edges::Vector)
  return DigraphWrap.IsMaximalMatching(GapObj(D), GapObj(edges, recursive = true))::Bool
end

@doc raw"""
    is_maximum_matching(D::Digraph, edges::Vector) -> Bool

Return `true` if `edges` is a maximum matching of the digraph `D`, and
`false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_maximum_matching(d, [[1, 2]])
true
```
"""
function is_maximum_matching(D::Digraph, edges::Vector)
  return DigraphWrap.IsMaximumMatching(GapObj(D), GapObj(edges, recursive = true))::Bool
end

@doc raw"""
    is_perfect_matching(D::Digraph, edges::Vector) -> Bool

Return `true` if `edges` is a perfect matching of the digraph `D`, i.e. a
matching in which every vertex is incident to an edge, and `false` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> is_perfect_matching(d, [[1, 2]])
true
```
"""
function is_perfect_matching(D::Digraph, edges::Vector)
  return DigraphWrap.IsPerfectMatching(GapObj(D), GapObj(edges, recursive = true))::Bool
end

@doc raw"""
    digraph_maximal_matching(D::Digraph) -> Vector{Tuple{Int,Int}}

Return a maximal matching of the digraph D.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_maximal_matching(d)
1-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
```
"""
function digraph_maximal_matching(D::Digraph)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphMaximalMatching(GapObj(D)))
end

@doc raw"""
    digraph_maximum_matching(D::Digraph) -> Vector{Tuple{Int,Int}}

Return a maximum matching of the digraph D.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_maximum_matching(d)
1-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
```
"""
function digraph_maximum_matching(D::Digraph)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphMaximumMatching(GapObj(D)))
end

@doc raw"""
    adjacency_matrix(D::Digraph) -> Matrix{Int}

Return the adjacency matrix of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> adjacency_matrix(d)
2×2 Matrix{Int64}:
 0  1
 1  0
```
"""
function adjacency_matrix(D::Digraph)
  return Matrix{Int}(DigraphWrap.AdjacencyMatrix(GapObj(D)))
end

@doc raw"""
    adjacency_matrix_mutable_copy(D::Digraph) -> Matrix{Int}

Return a mutable copy of the adjacency matrix of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> adjacency_matrix_mutable_copy(d)
2×2 Matrix{Int64}:
 0  1
 1  0
```
"""
function adjacency_matrix_mutable_copy(D::Digraph)
  return Matrix{Int}(DigraphWrap.AdjacencyMatrixMutableCopy(GapObj(D)))
end

@doc raw"""
    characteristic_polynomial(D::Digraph) -> String

Return the characteristic polynomial of the adjacency matrix of `D`.

# Examples
```jldoctest
julia> d = complete_digraph(2)
Digraph with 2 vertices, 2 edges

julia> characteristic_polynomial(d)
"x_1^2-1"
```
"""
function characteristic_polynomial(D::Digraph)
  return String(GAPWrap.String(DigraphWrap.CharacteristicPolynomial(GapObj(D))))
end

@doc raw"""
    boolean_adjacency_matrix(D::Digraph) -> Matrix{Bool}

Return the boolean adjacency matrix of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> boolean_adjacency_matrix(d)
2×2 Matrix{Bool}:
 0  1
 1  0
```
"""
function boolean_adjacency_matrix(D::Digraph)
  return Matrix{Bool}(DigraphWrap.BooleanAdjacencyMatrix(GapObj(D)))
end

@doc raw"""
    boolean_adjacency_matrix_mutable_copy(D::Digraph) -> Matrix{Bool}

Return a mutable copy of the boolean adjacency matrix of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> boolean_adjacency_matrix_mutable_copy(d)
2×2 Matrix{Bool}:
 0  1
 1  0
```
"""
function boolean_adjacency_matrix_mutable_copy(D::Digraph)
  return Matrix{Bool}(DigraphWrap.BooleanAdjacencyMatrixMutableCopy(GapObj(D)))
end

@doc raw"""
    digraph_adjacency_function(D::Digraph) -> GapObj

Return a GAP function which takes two integer parameters `x` and `y` and
returns `true` if there exists an edge from vertex `x` to vertex `y` in `D`,
and `false` otherwise.

# Examples
```jldoctest
julia> f = digraph_adjacency_function(digraph([[2], [1]]));

julia> f(1, 2)
true
```
"""
function digraph_adjacency_function(D::Digraph)
  return DigraphWrap.DigraphAdjacencyFunction(GapObj(D))
end

@doc raw"""
    digraph_source(D::Digraph) -> Vector{Int}

Return the source of each edge of the digraph `D` in order. Position `i` in
the output gives the source of the `i`th edge of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> digraph_source(d)
2-element Vector{Int64}:
 1
 2
```
"""
function digraph_source(D::Digraph)
  return [e[1] for e in digraph_edges(D)]
end

@doc raw"""
    digraph_range(D::Digraph) -> Vector{Int}

Return the range of each edge of the digraph `D` in order. Position `i` in
the output gives the range of the `i`th edge of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> digraph_range(d)
2-element Vector{Int64}:
 2
 1
```
"""
function digraph_range(D::Digraph)
  return [e[2] for e in digraph_edges(D)]
end

@doc raw"""
    out_neighbours(D::Digraph) -> Vector{Vector{Int}}

Return the out-neighbours of the digraph `D` as a list of lists.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], [1]])
Digraph with 3 vertices, 4 edges

julia> out_neighbours(d)
3-element Vector{Vector{Int64}}:
 [2, 3]
 [3]
 [1]
```
"""
function out_neighbours(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.OutNeighbours(GapObj(D)))
end

@doc raw"""
    out_neighbors(D::Digraph) -> Vector{Vector{Int}}

American spelling alias for [`out_neighbours`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> out_neighbors(d)
2-element Vector{Vector{Int64}}:
 [2]
 [1]
```
"""
function out_neighbors(D::Digraph)
  return out_neighbours(D)
end

@doc raw"""
    out_neighbours_mutable_copy(D::Digraph) -> Vector{Vector{Int}}

Return a mutable copy of the out-neighbours of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> out_neighbours_mutable_copy(d)
2-element Vector{Vector{Int64}}:
 [2]
 [1]
```
"""
function out_neighbours_mutable_copy(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.OutNeighboursMutableCopy(GapObj(D)))
end

@doc raw"""
    out_neighbors_mutable_copy(D::Digraph) -> Vector{Vector{Int}}

American spelling alias for [`out_neighbours_mutable_copy`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> out_neighbors_mutable_copy(d)
2-element Vector{Vector{Int64}}:
 [2]
 [1]
```
"""
function out_neighbors_mutable_copy(D::Digraph)
  return out_neighbours_mutable_copy(D)
end

@doc raw"""
    in_neighbours(D::Digraph) -> Vector{Vector{Int}}

Return the in-neighbours of the digraph `D` as a list of lists.

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges

julia> in_neighbours(d)
3-element Vector{Vector{Int64}}:
 [2]
 [1, 3]
 [2]
```
"""
function in_neighbours(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.InNeighbours(GapObj(D)))
end

@doc raw"""
    in_neighbors(D::Digraph) -> Vector{Vector{Int}}

American spelling alias for [`in_neighbours`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_neighbors(d)
2-element Vector{Vector{Int64}}:
 [2]
 [1]
```
"""
function in_neighbors(D::Digraph)
  return in_neighbours(D)
end

@doc raw"""
    in_neighbours_mutable_copy(D::Digraph) -> Vector{Vector{Int}}

Return a mutable copy of the in-neighbours of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_neighbours_mutable_copy(d)
2-element Vector{Vector{Int64}}:
 [2]
 [1]
```
"""
function in_neighbours_mutable_copy(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.InNeighboursMutableCopy(GapObj(D)))
end

@doc raw"""
    in_neighbors_mutable_copy(D::Digraph) -> Vector{Vector{Int}}

American spelling alias for [`in_neighbours_mutable_copy`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_neighbors_mutable_copy(d)
2-element Vector{Vector{Int64}}:
 [2]
 [1]
```
"""
function in_neighbors_mutable_copy(D::Digraph)
  return in_neighbours_mutable_copy(D)
end

@doc raw"""
    out_degrees(D::Digraph) -> Vector{Int}

Return the out-degrees of all vertices of `D`.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], [1]])
Digraph with 3 vertices, 4 edges

julia> out_degrees(d)
3-element Vector{Int64}:
 2
 1
 1
```
"""
function out_degrees(D::Digraph)
  return Vector{Int}(DigraphWrap.OutDegrees(GapObj(D)))
end

@doc raw"""
    out_degree_sequence(D::Digraph) -> Vector{Int}

Return the out-degrees of the vertices of the digraph `D`, sorted into
non-increasing order.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> out_degree_sequence(d)
2-element Vector{Int64}:
 1
 1
```
"""
function out_degree_sequence(D::Digraph)
  return Vector{Int}(DigraphWrap.OutDegreeSequence(GapObj(D)))
end

@doc raw"""
    out_degree_set(D::Digraph) -> Vector{Int}

Return the set of out-degrees of the vertices of the digraph `D`, sorted
into increasing order with duplicate entries removed.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> out_degree_set(d)
3-element Vector{Int64}:
 0
 1
 2
```
"""
function out_degree_set(D::Digraph)
  return Vector{Int}(DigraphWrap.OutDegreeSet(GapObj(D)))
end

@doc raw"""
    in_degrees(D::Digraph) -> Vector{Int}

Return the in-degrees of all vertices of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges

julia> in_degrees(d)
3-element Vector{Int64}:
 1
 2
 1
```
"""
function in_degrees(D::Digraph)
  return Vector{Int}(DigraphWrap.InDegrees(GapObj(D)))
end

@doc raw"""
    in_degree_sequence(D::Digraph) -> Vector{Int}

Return the in-degrees of the vertices of the digraph `D`, sorted into
non-increasing order.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_degree_sequence(d)
2-element Vector{Int64}:
 1
 1
```
"""
function in_degree_sequence(D::Digraph)
  return Vector{Int}(DigraphWrap.InDegreeSequence(GapObj(D)))
end

@doc raw"""
    in_degree_set(D::Digraph) -> Vector{Int}

Return the set of in-degrees of the vertices of the digraph `D`, sorted into
increasing order with duplicate entries removed.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> in_degree_set(d)
3-element Vector{Int64}:
 0
 1
 2
```
"""
function in_degree_set(D::Digraph)
  return Vector{Int}(DigraphWrap.InDegreeSet(GapObj(D)))
end

@doc raw"""
    out_degree_of_vertex(D::Digraph, v::Int) -> Int

Return the out-degree of vertex `v` in the digraph `D`, i.e. the number of
edges of `D` whose source is `v`.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> out_degree_of_vertex(d, 1)
2
```
"""
function out_degree_of_vertex(D::Digraph, v::Int)
  return Int(DigraphWrap.OutDegreeOfVertex(GapObj(D), v)::GapInt)
end

@doc raw"""
    out_neighbours_of_vertex(D::Digraph, v::Int) -> Vector{Int}

Return the out-neighbours of vertex `v` in the digraph `D`; a vertex may
appear more than once if there are multiple edges.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> out_neighbours_of_vertex(d, 1)
2-element Vector{Int64}:
 2
 2
```
"""
function out_neighbours_of_vertex(D::Digraph, v::Int)
  return Vector{Int}(DigraphWrap.OutNeighboursOfVertex(GapObj(D), v))
end

@doc raw"""
    out_neighbors_of_vertex(D::Digraph, v::Int) -> Vector{Int}

American spelling alias for [`out_neighbours_of_vertex`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1]])
Digraph with 2 vertices, 3 edges

julia> out_neighbors_of_vertex(d, 1)
2-element Vector{Int64}:
 2
 2
```
"""
function out_neighbors_of_vertex(D::Digraph, v::Int)
  return out_neighbours_of_vertex(D, v)
end

@doc raw"""
    in_degree_of_vertex(D::Digraph, v::Int) -> Int

Return the in-degree of vertex `v` in the digraph `D`, i.e. the number of
edges of `D` whose range is `v`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_degree_of_vertex(d, 1)
1
```
"""
function in_degree_of_vertex(D::Digraph, v::Int)
  return Int(DigraphWrap.InDegreeOfVertex(GapObj(D), v)::GapInt)
end

@doc raw"""
    in_neighbours_of_vertex(D::Digraph, v::Int) -> Vector{Int}

Return the in-neighbours of vertex `v` in the digraph `D`; a vertex may
appear more than once if there are multiple edges.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_neighbours_of_vertex(d, 1)
1-element Vector{Int64}:
 2
```
"""
function in_neighbours_of_vertex(D::Digraph, v::Int)
  return Vector{Int}(DigraphWrap.InNeighboursOfVertex(GapObj(D), v))
end

@doc raw"""
    in_neighbors_of_vertex(D::Digraph, v::Int) -> Vector{Int}

American spelling alias for [`in_neighbours_of_vertex`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> in_neighbors_of_vertex(d, 1)
1-element Vector{Int64}:
 2
```
"""
function in_neighbors_of_vertex(D::Digraph, v::Int)
  return in_neighbours_of_vertex(D, v)
end

@doc raw"""
    digraph_loops(D::Digraph) -> Vector{Int}

Return the list of the vertices of the digraph `D` at which there is a loop.

# Examples
```jldoctest
julia> d = digraph([[1, 2], [2]])
Digraph with 2 vertices, 3 edges

julia> digraph_loops(d)
2-element Vector{Int64}:
 1
 2
```
"""
function digraph_loops(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphLoops(GapObj(D)))
end

@doc raw"""
    degree_matrix(D::Digraph) -> Matrix{Int}

Return the degree matrix of `D` (diagonal matrix with out-degrees).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> degree_matrix(d)
2×2 Matrix{Int64}:
 1  0
 0  1
```
"""
function degree_matrix(D::Digraph)
  return Matrix{Int}(DigraphWrap.DegreeMatrix(GapObj(D)))
end

@doc raw"""
    laplacian_matrix(D::Digraph) -> Matrix{Int}

Return the Laplacian matrix of `D` (degree matrix minus adjacency matrix).

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> laplacian_matrix(d)
2×2 Matrix{Int64}:
  1  -1
 -1   1
```
"""
function laplacian_matrix(D::Digraph)
  return Matrix{Int}(DigraphWrap.LaplacianMatrix(GapObj(D)))
end

@doc raw"""
    partial_order_digraph_meet_of_vertices(D::Digraph, u::Int, v::Int) -> Int

If `D` is a partial order digraph, return the meet of the vertices `u` and
`v` of `D`. Throw an error if the meet does not exist.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> partial_order_digraph_meet_of_vertices(d, 2, 3)
2
```
"""
function partial_order_digraph_meet_of_vertices(D::Digraph, u::Int, v::Int)
  result = DigraphWrap.PartialOrderDigraphMeetOfVertices(GapObj(D), u, v)
  @req result != GAP.Globals.fail "partial_order_digraph_meet_of_vertices: the meet of u and v does not exist in D"
  return Int(result::GapInt)
end

@doc raw"""
    partial_order_digraph_join_of_vertices(D::Digraph, u::Int, v::Int) -> Int

If `D` is a partial order digraph, return the join of the vertices `u` and
`v` of `D`. Throw an error if the join does not exist.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 3], [3]])
Digraph with 3 vertices, 6 edges

julia> partial_order_digraph_join_of_vertices(d, 2, 3)
3
```
"""
function partial_order_digraph_join_of_vertices(D::Digraph, u::Int, v::Int)
  result = DigraphWrap.PartialOrderDigraphJoinOfVertices(GapObj(D), u, v)
  @req result != GAP.Globals.fail "partial_order_digraph_join_of_vertices: the join of u and v does not exist in D"
  return Int(result::GapInt)
end

@doc raw"""
    non_upper_semimodular_pair(D::Digraph) -> Tuple{Int, Int}

Return a pair of vertices of `D` witnessing the fact that `D` does not
represent an upper semimodular lattice. Throw an error if no such pair
exists (meaning that `D` is an upper semimodular lattice).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4, 5], [2, 4, 5], [3, 5], [4, 5], [5]])
Digraph with 5 vertices, 13 edges

julia> non_upper_semimodular_pair(d)
(3, 2)
```
"""
function non_upper_semimodular_pair(D::Digraph)
  result = DigraphWrap.NonUpperSemimodularPair(GapObj(D))
  @req result != GAP.Globals.fail "non_upper_semimodular_pair: no such pair exists (meaning that D is an upper semimodular lattice)"
  return Tuple{Int,Int}(result)
end
@doc raw"""
    non_lower_semimodular_pair(D::Digraph) -> Tuple{Int, Int}

Return a pair of vertices of `D` witnessing the fact that `D` does not
represent a lower semimodular lattice. Throw an error if no such pair
exists (meaning that `D` is a lower semimodular lattice).

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4, 5], [2, 4, 5], [3, 5], [4, 5], [5]])
Digraph with 5 vertices, 13 edges

julia> non_lower_semimodular_pair(d)
(4, 3)
```
"""
function non_lower_semimodular_pair(D::Digraph)
  result = DigraphWrap.NonLowerSemimodularPair(GapObj(D))
  @req result != GAP.Globals.fail "non_lower_semimodular_pair: no such pair exists (meaning that D is a lower semimodular lattice)"
  return Tuple{Int,Int}(result)
end

@doc raw"""
    digraph_diameter(D::Digraph) -> Int

Return the diameter of the digraph `D` (the greatest finite shortest
path distance between any pair of vertices).

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1, 2]])
Digraph with 3 vertices, 4 edges

julia> digraph_diameter(d)
2
```
"""
function digraph_diameter(D::Digraph)
  return Int(DigraphWrap.DigraphDiameter(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_shortest_distance(D::Digraph, s::Int, t::Int) -> Int

Return the shortest directed distance from `s` to `t` in `D`.
Throw an error if no such path exists.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_shortest_distance(d, 1, 3)
2
```
"""
function digraph_shortest_distance(D::Digraph, s::Int, t::Int)
  result = DigraphWrap.DigraphShortestDistance(GapObj(D), s, t)
  @req result != GAP.Globals.fail "digraph_shortest_distance: there is no directed path from s to t in D"
  return Int(result::GapInt)
end

@doc raw"""
    digraph_shortest_distance(D::Digraph, list::Vector{<:Integer}) -> Int

Return the length of the shortest directed path between the vertices
`list[1]` and `list[2]` of the digraph `D`. Throw an error if no such
path exists.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> digraph_shortest_distance(d, [1, 4])
2
```
"""
function digraph_shortest_distance(D::Digraph, list::Vector{<:Integer})
  result = DigraphWrap.DigraphShortestDistance(GapObj(D), list)
  @req result != GAP.Globals.fail "digraph_shortest_distance: there is no directed path between list[1] and list[2] in D"
  return Int(result::GapInt)
end

@doc raw"""
    digraph_shortest_distance(D::Digraph, list1::Vector{<:Integer}, list2::Vector{<:Integer}) -> Int

Return the length of the shortest directed path which starts at a vertex in
`list1` and terminates at a vertex in `list2`. Throw an error if no such
path exists; the result is `0` if `list1` and `list2` have non-empty
intersection.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> digraph_shortest_distance(d, [1, 2], [4])
1
```
"""
function digraph_shortest_distance(D::Digraph, list1::Vector{<:Integer}, list2::Vector{<:Integer})
  result = DigraphWrap.DigraphShortestDistance(GapObj(D), list1, list2)
  @req result != GAP.Globals.fail "digraph_shortest_distance: there is no directed path from a vertex in list1 to a vertex in list2 in D"
  return Int(result::GapInt)
end

@doc raw"""
    digraph_shortest_distances(D::Digraph) -> Matrix{Int}

Return the matrix of shortest directed distances between all pairs of
vertices of `D`. A value of `-1` indicates that no path exists.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1, 2]])
Digraph with 3 vertices, 4 edges

julia> shortest_distances(d)
3×3 Matrix{Int64}:
 0  1  2
 2  0  1
 1  1  0
```
"""
function shortest_distances(D::Digraph)
  result = DigraphWrap.DigraphShortestDistances(GapObj(D))
  n = nv(D)
  M = Matrix{Int}(undef, n, n)
  for i in 1:n
    for j in 1:n
      val = result[i][j]
      if val == GAP.Globals.fail
        M[i, j] = -1
      else
        M[i, j] = Int(val::GapInt)
      end
    end
  end
  return M
end

@doc raw"""
    digraph_longest_distance_from_vertex(D::Digraph, v::Int) -> Real

Return the length of the longest directed walk in the digraph `D` which
begins at vertex `v`. Return `Inf` if there is a directed walk of infinite
length starting at `v` (realised by repeatedly traversing a loop or directed
cycle).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_longest_distance_from_vertex(d, 1)
2
```
"""
function digraph_longest_distance_from_vertex(D::Digraph, v::Int)
  result = DigraphWrap.DigraphLongestDistanceFromVertex(GapObj(D), v)
  result == GAP.Globals.infinity && return Inf
  return Int(result::GapInt)
end

@doc raw"""
    digraph_distance_set(D::Digraph, v::Int, distance::Int) -> Vector{Int}

    digraph_distance_set(D::Digraph, v::Int, distances::Vector{<:Integer}) -> Vector{Int}

Return the list of all vertices of the digraph `D` whose shortest distance to
vertex `v` is equal to `distance`, or is contained in the list `distances`.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> digraph_distance_set(d, 1, 2)
1-element Vector{Int64}:
 4

julia> digraph_distance_set(d, 1, [1, 2])
3-element Vector{Int64}:
 2
 3
 4
```
"""
function digraph_distance_set(D::Digraph, v::Int, distance::Int)
  return Vector{Int}(DigraphWrap.DigraphDistanceSet(GapObj(D), v, distance))
end

function digraph_distance_set(D::Digraph, v::Int, distances::Vector{<:Integer})
  return Vector{Int}(DigraphWrap.DigraphDistanceSet(GapObj(D), v, distances))
end

@doc raw"""
    digraph_girth(D::Digraph) -> Real

Return the girth (the length of the shortest directed cycle) of the digraph
`D`. Return `Inf` if `D` has no directed cycles.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_girth(d)
3
```
"""
function digraph_girth(D::Digraph)
  result = DigraphWrap.DigraphGirth(GapObj(D))
  result == GAP.Globals.infinity && return Inf
  return Int(result::GapInt)
end

@doc raw"""
    digraph_odd_girth(D::Digraph) -> Real

Return the odd girth (the length of the shortest directed cycle of odd
length) of the digraph `D`. Return `Inf` if `D` has no directed cycles of
odd length.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_odd_girth(d)
3
```
"""
function digraph_odd_girth(D::Digraph)
  result = DigraphWrap.DigraphOddGirth(GapObj(D))
  result == GAP.Globals.infinity && return Inf
  return Int(result::GapInt)
end

@doc raw"""
    digraph_undirected_girth(D::Digraph) -> Real

If `D` is a symmetric digraph, return the girth of `D` when treated as an
undirected graph (each pair of edges `[i, j]` and `[j, i]` is treated as a
single edge). Return `Inf` if `D` has no undirected cycles.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> digraph_undirected_girth(d)
3
```
"""
function digraph_undirected_girth(D::Digraph)
  result = DigraphWrap.DigraphUndirectedGirth(GapObj(D))
  result == GAP.Globals.infinity && return Inf
  return Int(result::GapInt)
end

@doc raw"""
    digraph_connected_components(D::Digraph) -> Dict

Return the connected components of `D` as a GAP record with components
`comps` (list of vertex sets) and `id` (component id of each vertex).

# Examples
```jldoctest
julia> d = digraph([[2], [1], []])
Digraph with 3 vertices, 2 edges

julia> digraph_connected_components(d)
Dict{Symbol, Any} with 2 entries:
  :id    => Any[1, 1, 2]
  :comps => Any[Any[1, 2], Any[3]]
```
"""
function digraph_connected_components(D::Digraph)
  d = DigraphWrap.DigraphConnectedComponents(GapObj(D))
  return Dict{Symbol, Any}(d, recursive = true)
end

@doc raw"""
    digraph_nr_connected_components(D::Digraph) -> Int

Return the number of connected components of the digraph D.

# Examples
```jldoctest
julia> d = digraph([[2], [1], [4], [3]])
Digraph with 4 vertices, 4 edges

julia> digraph_nr_connected_components(d)
2
```
"""
function digraph_nr_connected_components(D::Digraph)
  return Int(DigraphWrap.DigraphNrConnectedComponents(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_connected_component(D::Digraph, v::Int) -> Vector{Int}

Return the vertices in the connected component of `D` containing vertex `v`.

# Examples
```jldoctest
julia> d = digraph([[2], [1], []])
Digraph with 3 vertices, 2 edges

julia> digraph_connected_component(d, 1)
2-element Vector{Int64}:
 1
 2
```
"""
function digraph_connected_component(D::Digraph, v::Int)
  return Vector{Int}(DigraphWrap.DigraphConnectedComponent(GapObj(D), v))
end

@doc raw"""
    digraph_strongly_connected_components(D::Digraph) -> Dict

Return the strongly connected components of `D` as a GAP record.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1], []])
Digraph with 4 vertices, 3 edges

julia> digraph_strongly_connected_components(d)
Dict{Symbol, Any} with 2 entries:
  :id    => Any[1, 1, 1, 2]
  :comps => Any[Any[1, 2, 3], Any[4]]
```
"""
function digraph_strongly_connected_components(D::Digraph)
  d = DigraphWrap.DigraphStronglyConnectedComponents(GapObj(D))
  return Dict{Symbol, Any}(d, recursive = true)
end

@doc raw"""
    digraph_nr_strongly_connected_components(D::Digraph) -> Int

Return the number of strongly connected components of the digraph D.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_nr_strongly_connected_components(d)
3
```
"""
function digraph_nr_strongly_connected_components(D::Digraph)
  return Int(DigraphWrap.DigraphNrStronglyConnectedComponents(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_strongly_connected_component(D::Digraph, v::Int) -> Vector{Int}

Return the vertices in the strongly connected component of `D` containing
vertex `v`.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_strongly_connected_component(d, 1)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function digraph_strongly_connected_component(D::Digraph, v::Int)
  return Vector{Int}(DigraphWrap.DigraphStronglyConnectedComponent(GapObj(D), v))
end

@doc raw"""
    digraph_bicomponents(D::Digraph) -> Tuple{Vector{Int}, Vector{Int}}

Return the bicomponents (bipartition) of the digraph `D`. Throw an error if `D` is not bipartite.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1], [1]])
Digraph with 3 vertices, 4 edges

julia> digraph_bicomponents(d)
([1], [2, 3])
```
"""
function digraph_bicomponents(D::Digraph)
  result = DigraphWrap.DigraphBicomponents(GapObj(D))
  @req result != GAP.Globals.fail "digraph_bicomponents: the digraph is not bipartite"
  return Tuple{Vector{Int}, Vector{Int}}(result)
end

@doc raw"""
    articulation_points(D::Digraph) -> Vector{Int}

Return the list of articulation points (cut vertices) of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges

julia> articulation_points(d)
1-element Vector{Int64}:
 2
```
"""
function articulation_points(D::Digraph)
  return Vector{Int}(DigraphWrap.ArticulationPoints(GapObj(D)))
end

@doc raw"""
    minimal_cyclic_edge_cut(D::Digraph) -> Vector{Int}

Return a minimal set of edges whose removal makes D acyclic.

# Examples
```jldoctest
julia> d = hypercube_graph(3)
Digraph with 8 vertices, 24 edges

julia> minimal_cyclic_edge_cut(d)
4-element Vector{Vector{Int64}}:
 [1, 5]
 [2, 6]
 [4, 8]
 [3, 7]
```
"""
function minimal_cyclic_edge_cut(D::Digraph)
  result = DigraphWrap.MinimalCyclicEdgeCut(GapObj(D))
  @req result != GAP.Globals.fail "The digraph is not cubic, not connected, has less than 8 vertices or does not have a cyclic edge cut."
  return Vector{Vector{Int64}}(result)
end

@doc raw"""
    bridges(D::Digraph) -> Vector{Tuple{Int,Int}}

Return the list of bridges of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges

julia> bridges(d)
2-element Vector{Tuple{Int64, Int64}}:
 (2, 3)
 (1, 2)
```
"""
function bridges(D::Digraph)
  return Vector{Tuple{Int,Int}}(DigraphWrap.Bridges(GapObj(D)))
end

@doc raw"""
    strong_orientation(D::Digraph) -> Digraph

Return a strong orientation of the digraph D, if D is bridgeless.

# Examples
```jldoctest
julia> d = cycle_symmetric_digraph(3)
Digraph with 3 vertices, 6 edges

julia> strong_orientation(d)
Digraph with 3 vertices, 3 edges
```
"""
function strong_orientation(D::Digraph)
  result = DigraphWrap.StrongOrientation(GapObj(D))
  @req result != GAP.Globals.fail "D is symmetric but does not admit a strong orientation."
  return Digraph(result)
end

@doc raw"""
    strong_orientation_attr(D::Digraph) -> Digraph

Return a strong orientation of the digraph `D` as an attribute value.
Throw an error if `D` does not admit a strong orientation.

# Examples
```jldoctest
julia> d = cycle_symmetric_digraph(3)
Digraph with 3 vertices, 6 edges

julia> strong_orientation_attr(d) isa Digraph
true
```
"""
function strong_orientation_attr(D::Digraph)
  result = DigraphWrap.StrongOrientationAttr(GapObj(D))
  @req result != GAP.Globals.fail "strong_orientation_attr: the digraph does not admit a strong orientation"
  return Digraph(result)
end

@doc raw"""
    digraph_period(D::Digraph) -> Int

Return the period (gcd of directed cycle lengths) of the digraph `D`.
Returns `0` if the digraph has no directed cycles.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_period(d)
3
```
"""
function digraph_period(D::Digraph)
  return Int(DigraphWrap.DigraphPeriod(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_floyd_warshall(D::Digraph, func, nopath, edge) -> GapObj

Run a generalised version of the Floyd-Warshall algorithm on the digraph `D`
and return the resulting matrix as a GAP object. The matrix is initialised
by setting entry `[i][j]` to `edge` if there is an edge with source `i` and
range `j`, and to `nopath` otherwise; the function `func` is then called with
the matrix and three positive integers inside three nested loops over the
vertices of `D`.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> digraph_floyd_warshall(d, (mat, i, j, k) -> nothing, 0, 1) isa GapObj
true
```
"""
function digraph_floyd_warshall(D::Digraph, func, nopath, edge)
  func isa Function && (func = _gap_wrap_function(func))
  return DigraphWrap.DigraphFloydWarshall(GapObj(D), func, nopath, edge)
end

const _GAP_FUNC_WRAPPER = Ref{Any}(nothing)

# Convert a Julia function to a GAP function by wrapping it in a GAP lambda,
# so that GAP operations whose methods filter on `IsFunction` accept it.
function _gap_wrap_function(func)
  if _GAP_FUNC_WRAPPER[] === nothing
    _GAP_FUNC_WRAPPER[] = GAP.evalstr("f -> function(arg) return CallFuncList(f, arg); end")
  end
  return _GAP_FUNC_WRAPPER[](func)
end

@doc raw"""
    is_reachable(D::Digraph, s::Int, t::Int) -> Bool

Return `true` if there is a directed path from `s` to `t` in `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> is_reachable(d, 1, 3)
true

julia> is_reachable(d, 3, 1)
false
```
"""
function is_reachable(D::Digraph, s::Int, t::Int)
  return DigraphWrap.IsReachable(GapObj(D), s, t)::Bool
end

@doc raw"""
    is_digraph_path(D::Digraph, v::Vector{<:Integer}, a::Vector{<:Integer}) -> Bool

    is_digraph_path(D::Digraph, list::Vector{<:Vector{<:Integer}}) -> Bool

Return `true` if `v` and `a` describe a directed path (or directed cycle) in
the digraph `D`: `v` is the list of vertices of the path and `a` is the list
of positions of each successive vertex in the out-neighbours of the previous
vertex. The second form is equivalent to calling the first form with
`list[1]` and `list[2]`.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> is_digraph_path(d, [1, 2, 4], [1, 1])
true

julia> is_digraph_path(d, [[1, 2, 4], [1, 1]])
true
```
"""
function is_digraph_path(D::Digraph, v::Vector{<:Integer}, a::Vector{<:Integer})
  return DigraphWrap.IsDigraphPath(GapObj(D), GapObj(v, recursive = true), GapObj(a, recursive = true))::Bool
end

function is_digraph_path(D::Digraph, list::Vector{<:Vector{<:Integer}})
  return DigraphWrap.IsDigraphPath(GapObj(D), GapObj(list, recursive = true))::Bool
end

@doc raw"""
    vertices_reachable_from(D::Digraph, root::Int) -> Vector{Int}

    vertices_reachable_from(D::Digraph, list::Vector{<:Integer}) -> Vector{Int}

Return a list of the vertices `v` for which there exists a non-trivial
directed walk from `root` (or from any vertex in `list`) to `v` in the
digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> vertices_reachable_from(d, 1)
3-element Vector{Int64}:
 2
 3
 4
```
"""
function vertices_reachable_from(D::Digraph, root::Int)
  return Vector{Int}(DigraphWrap.VerticesReachableFrom(GapObj(D), root))
end

function vertices_reachable_from(D::Digraph, list::Vector{<:Integer})
  return Vector{Int}(DigraphWrap.VerticesReachableFrom(GapObj(D), list))
end

@doc raw"""
    digraph_path(D::Digraph, s::Int, t::Int) -> Vector{Vector{Int}}

Return a directed path from `s` to `t` in `D` as a list
`[vertices, edges]` where `vertices` lists the vertices along the
path and `edges` contains edge indices (positions in
`OutNeighbours(D)`). Throw an error if no path exists.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [2, 3]])
Digraph with 3 vertices, 4 edges

julia> digraph_path(d, 1, 3)
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 1]
```
"""
function digraph_path(D::Digraph, s::Int, t::Int)
  result = DigraphWrap.DigraphPath(GapObj(D), s, t)
  @req result != GAP.Globals.fail "digraph_path: there is no directed path from s to t in D"
  return Vector{Vector{Int64}}(result)
end

@doc raw"""
    digraph_shortest_path(D::Digraph, s::Int, t::Int) -> Vector{Vector{Int}}

Return a shortest directed path from `s` to `t` in `D` as a list
`[vertices, edges]` (same format as `digraph_path`).
Throw an error if no path exists.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [2, 3]])
Digraph with 3 vertices, 4 edges

julia> digraph_shortest_path(d, 1, 3)
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 1]
```
"""
function digraph_shortest_path(D::Digraph, s::Int, t::Int)
  result = DigraphWrap.DigraphShortestPath(GapObj(D), s, t)
  @req result != GAP.Globals.fail "digraph_shortest_path: there is no directed path from s to t in D"
  return Vector{Vector{Int64}}(result)
end

@doc raw"""
    digraph_random_walk(D::Digraph, v::Int, t::Int) -> Vector{Vector{Int}}

Return a directed path corresponding to a random walk in the digraph `D`,
starting at vertex `v` and having length no more than `t`. The output has
the same form as the output of [`digraph_path`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_random_walk(d, 1, 5)
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 1]
```
"""
function digraph_random_walk(D::Digraph, v::Int, t::Int)
  return Vector{Vector{Int64}}(DigraphWrap.DigraphRandomWalk(GapObj(D), v, t))
end

@doc raw"""
    digraph_absorption_probabilities(D::Digraph) -> Vector{Vector{Float64}}

Return the absorption probabilities of the absorbing Markov chain of D.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_absorption_probabilities(d)
3-element Vector{Vector{Float64}}:
 [1.0, 0.0, 0.0]
 [1.0, 0.0, 0.0]
 [1.0, 0.0, 0.0]
```
"""
function digraph_absorption_probabilities(D::Digraph)
  return Vector{Vector{Float64}}(DigraphWrap.DigraphAbsorptionProbabilities(GapObj(D)))
end

@doc raw"""
    digraph_absorption_expected_steps(D::Digraph) -> Vector{Float64}

Return the expected number of steps before absorption in the absorbing
Markov chain of D.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_absorption_expected_steps(d)
3-element Vector{Float64}:
 2.0
 1.0
 0.0
```
"""
function digraph_absorption_expected_steps(D::Digraph)
  return Vector{Float64}(DigraphWrap.DigraphAbsorptionExpectedSteps(GapObj(D)))
end

@doc raw"""
    dominators(D::Digraph, root::Int) -> Vector{Vector{Int}}

Return the dominators of each vertex of the digraph `D` with respect to
`root`: position `i` of the output contains the vertices that are contained
in every path from `root` to vertex `i`, not including `i` itself. If there
is no path from `root` to vertex `i`, then position `i` is the empty list.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> dominators(d, 1)
4-element Vector{Vector{Int64}}:
 []
 [1]
 [1]
 [1]
```
"""
function dominators(D::Digraph, root::Int)
  result = DigraphWrap.Dominators(GapObj(D), root)
  n = length(result)
  out = Vector{Vector{Int}}(undef, n)
  for i in 1:n
    if GAP.Wrappers.ISB_LIST(result, i)
      out[i] = Vector{Int}(result[i])
    else
      out[i] = Int[]
    end
  end
  return out
end

@doc raw"""
    dominator_tree(D::Digraph, root::Int) -> NamedTuple

Return the dominator tree of the digraph `D` with respect to `root` as a
named tuple with fields `:idom` (the immediate dominators of the vertices)
and `:preorder` (the preorder values of the vertices defined by the depth
first search executed on `D`).

# Examples
```jldoctest
julia> d = digraph([[2, 3], [4], [4], []])
Digraph with 4 vertices, 4 edges

julia> dominator_tree(d, 1)
(idom = Any[nothing, 1, 1, 1], preorder = [1, 2, 4, 3])
```
"""
function dominator_tree(D::Digraph, root::Int)
  result = DigraphWrap.DominatorTree(GapObj(D), root)
  idom = Any[i == GAP.Globals.fail ? nothing : Int(i::GapInt) for i in result.idom]
  preorder = Vector{Int}(result.preorder)
  return (idom = idom, preorder = preorder)
end
@doc raw"""
    iterator_of_paths(D::Digraph, u::Int, v::Int) -> GapObj

Return a GAP iterator of the non-trivial directed paths (or directed cycles,
in the case `u == v`) in the digraph `D` from vertex `u` to vertex `v`. Use
`GAP.Globals.NextIterator` to obtain the next path and
`GAP.Globals.IsDoneIterator` to test whether all paths have been returned.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> iterator_of_paths(d, 1, 2) isa GapObj
true
```
"""
function iterator_of_paths(D::Digraph, u::Int, v::Int)
  return DigraphWrap.IteratorOfPaths(GapObj(D), u, v)
end

@doc raw"""
    digraph_all_simple_circuits(D::Digraph) -> Vector{Vector{Int}}

Return a list of all simple directed circuits of `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_all_simple_circuits(d)
1-element Vector{Vector{Int64}}:
 [1, 2, 3]
```
"""
function digraph_all_simple_circuits(D::Digraph)
  return Vector{Vector{Int64}}(DigraphWrap.DigraphAllSimpleCircuits(GapObj(D)))
end

@doc raw"""
    digraph_longest_simple_circuit(D::Digraph) -> Vector{Int}

Return a longest simple circuit of the digraph `D`. Throw an error if `D`
has no simple circuits.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_longest_simple_circuit(d)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function digraph_longest_simple_circuit(D::Digraph)
  result = DigraphWrap.DigraphLongestSimpleCircuit(GapObj(D))
  @req result != GAP.Globals.fail "digraph_longest_simple_circuit: the digraph has no simple circuits"
  return Vector{Int}(result)
end

@doc raw"""
    digraph_all_undirected_simple_circuits(D::Digraph) -> Vector{Vector{Int}}

Return a list of the undirected simple circuits of the digraph `D`. A simple
circuit is undirected if the orientation of the edges in the circuit does not
matter. Multiple edges are ignored.

# Examples
```jldoctest
julia> digraph_all_undirected_simple_circuits(cycle_digraph(4))
1-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
```
"""
function digraph_all_undirected_simple_circuits(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.DigraphAllUndirectedSimpleCircuits(GapObj(D)))
end

@doc raw"""
    digraph_all_chordless_cycles(D::Digraph) -> Vector{Vector{Int}}

Return a list of the chordless cycles of the digraph `D`, i.e. undirected
simple circuits in which no pair of vertices is connected by an edge that is
not in the circuit. Cycles of length two are ignored.

# Examples
```jldoctest
julia> digraph_all_chordless_cycles(cycle_digraph(4))
1-element Vector{Vector{Int64}}:
 [2, 1, 4, 3]
```
"""
function digraph_all_chordless_cycles(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.DigraphAllChordlessCycles(GapObj(D)))
end

@doc raw"""
    digraph_all_chordless_cycles_of_maximal_length(D::Digraph, max_length::Int) -> Vector{Vector{Int}}

Return a list of all chordless cycles of the digraph `D` of length at most
`max_length`.

# Examples
```jldoctest
julia> digraph_all_chordless_cycles_of_maximal_length(cycle_digraph(4), 4)
1-element Vector{Vector{Int64}}:
 [2, 1, 4, 3]
```
"""
function digraph_all_chordless_cycles_of_maximal_length(D::Digraph, max_length::Int)
  return Vector{Vector{Int}}(DigraphWrap.DigraphAllChordlessCyclesOfMaximalLength(GapObj(D), max_length))
end

@doc raw"""
    facial_walks(D::Digraph, list::Vector{<:Vector{<:Integer}}) -> Vector{Vector{Int}}

If `D` is an Eulerian digraph and `list` is a rotation system of `D`, return
a list of the facial walks in `D`. Multiple edges and loops are ignored.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> facial_walks(d, [[2, 3], [1, 3], [1, 2]])
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 3, 2]
```
"""
function facial_walks(D::Digraph, list::Vector{<:Vector{<:Integer}})
  return Vector{Vector{Int}}(DigraphWrap.FacialWalks(GapObj(D), GapObj(list, recursive = true)))
end

@doc raw"""
    digraph_layers(D::Digraph, v::Int) -> Vector{Vector{Int}}

Return a list of layers of `D` starting from vertex `v`.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_layers(d, 1)
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [3]
```
"""
function digraph_layers(D::Digraph, v::Int)
  return Vector{Vector{Int}}(DigraphWrap.DigraphLayers(GapObj(D), v))
end

@doc raw"""
    digraph_degeneracy(D::Digraph) -> Int

Return the degeneracy of the digraph `D`. Throw an error if `D` is not
symmetric or has multiple edges.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_degeneracy(d)
2
```
"""
function digraph_degeneracy(D::Digraph)
  result = DigraphWrap.DigraphDegeneracy(GapObj(D))
  @req result != GAP.Globals.fail "digraph_degeneracy: the digraph is not symmetric or has multiple edges, so its degeneracy is not defined"
  return Int(result::GapInt)
end

@doc raw"""
    digraph_degeneracy_ordering(D::Digraph) -> Vector{Int}

Return a degeneracy ordering of the vertices of D.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_degeneracy_ordering(d)
3-element Vector{Int64}:
 3
 2
 1
```
"""
function digraph_degeneracy_ordering(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphDegeneracyOrdering(GapObj(D)))
end

@doc raw"""
    hamiltonian_path(D::Digraph) -> Vector{Int}

Return a Hamiltonian path in the digraph `D`. Throw an error if `D` has no
Hamiltonian path.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> hamiltonian_path(d)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function hamiltonian_path(D::Digraph)
  result = DigraphWrap.HamiltonianPath(GapObj(D))
  @req result != GAP.Globals.fail "hamiltonian_path: the digraph has no Hamiltonian path"
  return Vector{Int}(result)
end

@doc raw"""
    nr_spanning_trees(D::Digraph) -> Int

Return the number of spanning trees of the digraph `D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> nr_spanning_trees(d)
3
```
"""
function nr_spanning_trees(D::Digraph)
  return Int(DigraphWrap.NrSpanningTrees(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_dijkstra(D::Digraph, s::Int) -> Vector{Vector{Int}}

    digraph_dijkstra(D::Digraph, s::Int, t::Int) -> Vector{Vector{Int}}

Run Dijkstra's algorithm on the digraph `D` from source vertex `s`. Return a
pair of lists: the first list contains the distances of the vertices from
`s` (with `-1` for vertices that are not reachable from `s`), and the second
list contains the previous vertex in the shortest path from `s` to each
vertex (with `-1` for `s` and for unvisited vertices). If the optional
target vertex `t` is given, then the computation stops once the shortest
distance to `t` is known.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_dijkstra(d, 1)
2-element Vector{Vector{Int64}}:
 [0, 1, 2]
 [-1, 1, 2]

julia> digraph_dijkstra(d, 1, 3)
2-element Vector{Vector{Int64}}:
 [0, 1, 2]
 [-1, 1, 2]
```
"""
function digraph_dijkstra(D::Digraph, s::Int)
  return _dijkstra_result(DigraphWrap.DigraphDijkstra(GapObj(D), s))
end

function digraph_dijkstra(D::Digraph, s::Int, t::Int)
  return _dijkstra_result(DigraphWrap.DigraphDijkstra(GapObj(D), s, t))
end

# Convert the output of DigraphDijkstra (a pair of GAP lists) to a pair of
# Julia vectors of integers, mapping GAP `infinity` to `-1`.
function _dijkstra_result(result)
  dists = Int[]
  for i in 1:length(result[1])
    val = result[1][i]
    push!(dists, val == GAP.Globals.infinity ? -1 : Int(val::GapInt))
  end
  prev = Int[]
  for i in 1:length(result[2])
    push!(prev, Int(result[2][i]::GapInt))
  end
  return [dists, prev]
end

@doc raw"""
    digraph_vertex_connectivity(D::Digraph) -> Int

Return the vertex connectivity of the digraph `D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_vertex_connectivity(d)
2
```
"""
function digraph_vertex_connectivity(D::Digraph)
  return Int(DigraphWrap.DigraphVertexConnectivity(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_cycle_basis(D::Digraph) -> GapObj

Return a cycle basis of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_cycle_basis(d)
([[2], [3], [1]], GAP: [ <a GF2 vector of length 3> ])
```
"""
function digraph_cycle_basis(D::Digraph)
  result = DigraphWrap.DigraphCycleBasis(GapObj(D))
  return (Vector{Vector{Int64}}(result[1]), result[2])
end

@doc raw"""
    digraph_is_king(D::Digraph, v::Int, k::Int) -> Bool

If `D` is a tournament and `v` is a vertex of `D`, return `true` if every
other vertex of `D` is reachable from `v` by a path of length at most `k`,
and `false` otherwise.

# Examples
```jldoctest
julia> t = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> digraph_is_king(t, 1, 1)
true
```
"""
function digraph_is_king(D::Digraph, v::Int, k::Int)
  return DigraphWrap.DigraphIsKing(GapObj(D), v, k)::Bool
end

@doc raw"""
    digraph_kings(D::Digraph, n::Int) -> Vector{Int}

If `D` is a tournament, return a list of the `n`-kings of `D` (see
[`digraph_is_king`](@ref)).

# Examples
```jldoctest
julia> t = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> digraph_kings(t, 1)
1-element Vector{Int64}:
 1
```
"""
function digraph_kings(D::Digraph, n::Int)
  return Vector{Int}(DigraphWrap.DigraphKings(GapObj(D), n))
end

@doc raw"""
    group_of_cayley_digraph(D::Digraph) -> PermGroup

Return the group of which `D` is a Cayley digraph, if applicable.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3));

julia> group_of_cayley_digraph(d)
Symmetric group of degree 3
```
"""
function group_of_cayley_digraph(D::Digraph)
  return PermGroup(DigraphWrap.GroupOfCayleyDigraph(GapObj(D)))
end

@doc raw"""
    semigroup_of_cayley_digraph(D::Digraph) -> GapObj

If `D` is a Cayley digraph of a semigroup `S`, return `S`.

# Examples
```jldoctest
julia> D = cayley_digraph(symmetric_group(3))
Digraph with 6 vertices, 12 edges

julia> semigroup_of_cayley_digraph(D) isa GapObj
true
```
"""
function semigroup_of_cayley_digraph(D::Digraph)
  return DigraphWrap.SemigroupOfCayleyDigraph(GapObj(D))
end

@doc raw"""
    generators_of_cayley_digraph(D::Digraph) -> GapObj

Return the generating set of the Cayley digraph `D`, if applicable.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3));

julia> generators_of_cayley_digraph(d)
2-element Vector{PermGroupElem}:
 (1,2,3)
 (1,2)
```
"""
function generators_of_cayley_digraph(D::Digraph)
  return gens(group_of_cayley_digraph(D))
end

@doc raw"""
    as_semigroup(filt, D::Digraph) -> GapObj

If `filt` is equal to `GAP.Globals.IsPartialPermSemigroup` and `D` is a join
semilattice or lattice digraph, return a semigroup of partial permutations
which is isomorphic to the semigroup whose elements are the vertices of `D`
with the binary operation `partial_order_digraph_join_of_vertices`. If `D`
is a meet semilattice but not a join semilattice, then the meet operation is
used instead.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2], [3]])
Digraph with 3 vertices, 5 edges

julia> as_semigroup(GAP.Globals.IsPartialPermSemigroup, d) isa GapObj
true
```
"""
function as_semigroup(filt, D::Digraph)
  return DigraphWrap.AsSemigroup(filt, GapObj(D))
end

@doc raw"""
    as_monoid(filt, D::Digraph) -> GapObj

If `D` is a lattice digraph, return a monoid of partial permutations which is
isomorphic to the monoid whose elements are the vertices of `D` with the
binary operation `partial_order_digraph_join_of_vertices`.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 4], [2, 3, 4], [3, 4], [4]])
Digraph with 4 vertices, 10 edges

julia> as_monoid(GAP.Globals.IsPartialPermMonoid, d) isa GapObj
true
```
"""
function as_monoid(filt, D::Digraph)
  return DigraphWrap.AsMonoid(filt, GapObj(D))
end

@doc raw"""
    as_semigroup(filt, Y::Digraph, gps, homs) -> GapObj

If `filt` is equal to `GAP.Globals.IsPartialPermSemigroup`, `Y` is a join or
meet semilattice digraph, `gps` is a list of groups corresponding to each
vertex of `Y`, and `homs` is a list containing a triple `[i, j, hom]` for
each edge `(i, j)` in the transitive reduction of `Y`, where `hom` is a group
homomorphism from `gps[i]` to `gps[j]`, then return a semigroup of partial
permutations which is isomorphic to the strong semilattice of groups.

# Examples
```jldoctest
julia> Y = digraph([[1, 2], [2]])
Digraph with 2 vertices, 3 edges

julia> g = GAP.Globals.CyclicGroup(2);

julia> homs = GapObj([[1, 2, GAP.Globals.IdentityMapping(g)]], recursive=true);

julia> as_semigroup(GAP.Globals.IsPartialPermSemigroup, Y, [g, g], homs) isa GapObj
true
```
"""
function as_semigroup(filt, Y::Digraph, gps, homs)
  return DigraphWrap.AsSemigroup(filt, GapObj(Y), gps, homs)
end

@doc raw"""
    kuratowski_planar_subdigraph(D::Digraph) -> Vector{Vector{Int}}

Return the immutable list of lists of out-neighbours of an induced
subdigraph of `D` (excluding multiple edges and loops) that witnesses the
fact that `D` is not planar. Throw an error if `D` is planar.

# Examples
```jldoctest
julia> kuratowski_planar_subdigraph(complete_digraph(5)) isa Vector{Vector{Int}}
true
```
"""
function kuratowski_planar_subdigraph(D::Digraph)
  result = DigraphWrap.KuratowskiPlanarSubdigraph(GapObj(D))
  @req result != GAP.Globals.fail "kuratowski_planar_subdigraph: the digraph is planar, so no witness for non-planarity exists"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    kuratowski_outer_planar_subdigraph(D::Digraph) -> Vector{Vector{Int}}

Return the immutable list of lists of out-neighbours of an induced
subdigraph of `D` (excluding multiple edges and loops) that witnesses the
fact that `D` is not outer planar. Throw an error if `D` is outer planar.

# Examples
```jldoctest
julia> kuratowski_outer_planar_subdigraph(complete_digraph(4)) isa Vector{Vector{Int}}
true
```
"""
function kuratowski_outer_planar_subdigraph(D::Digraph)
  result = DigraphWrap.KuratowskiOuterPlanarSubdigraph(GapObj(D))
  @req result != GAP.Globals.fail "kuratowski_outer_planar_subdigraph: the digraph is outer planar, so no witness for non-outer-planarity exists"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    planar_embedding(D::Digraph) -> Vector{Vector{Int}}

If `D` is a planar digraph, return the immutable list of lists of
out-neighbours of `D` (excluding multiple edges and loops) such that each
vertex's neighbours are given in clockwise order. Throw an error if `D` is
not planar.

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2, 4], [1, 3]])
Digraph with 4 vertices, 7 edges

julia> planar_embedding(d)
4-element Vector{Vector{Int64}}:
 [2]
 [3, 1]
 [4, 2]
 [1, 3]
```
"""
function planar_embedding(D::Digraph)
  result = DigraphWrap.PlanarEmbedding(GapObj(D))
  @req result != GAP.Globals.fail "planar_embedding: the digraph is not planar"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    outer_planar_embedding(D::Digraph) -> Vector{Vector{Int}}

If `D` is an outer planar digraph, return the immutable list of lists of
out-neighbours of `D` (excluding multiple edges and loops) such that each
vertex's neighbours are given in clockwise order. Throw an error if `D` is
not outer planar.

# Examples
```jldoctest
julia> d = digraph([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges

julia> outer_planar_embedding(d)
3-element Vector{Vector{Int64}}:
 [2]
 [1, 3]
 [2]
```
"""
function outer_planar_embedding(D::Digraph)
  result = DigraphWrap.OuterPlanarEmbedding(GapObj(D))
  @req result != GAP.Globals.fail "outer_planar_embedding: the digraph is not outer planar"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    subdigraph_homeomorphic_to_k23(D::Digraph) -> Vector{Vector{Int}}

Return the immutable list of lists of out-neighbours of a subdigraph of `D`
which is homeomorphic to the complete bipartite graph with vertex sets of
sizes 2 and 3. Throw an error if `D` has no such subdigraph.

# Examples
```jldoctest
julia> subdigraph_homeomorphic_to_k23(complete_digraph(5)) isa Vector{Vector{Int}}
true
```
"""
function subdigraph_homeomorphic_to_k23(D::Digraph)
  result = DigraphWrap.SubdigraphHomeomorphicToK23(GapObj(D))
  @req result != GAP.Globals.fail "subdigraph_homeomorphic_to_k23: the digraph has no subdigraph homeomorphic to K23"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    subdigraph_homeomorphic_to_k33(D::Digraph) -> Vector{Vector{Int}}

Return the immutable list of lists of out-neighbours of a subdigraph of `D`
which is homeomorphic to the complete bipartite graph with vertex sets of
sizes 3 and 3. Throw an error if `D` has no such subdigraph.

# Examples
```jldoctest
julia> subdigraph_homeomorphic_to_k33(complete_digraph(6)) isa Vector{Vector{Int}}
true
```
"""
function subdigraph_homeomorphic_to_k33(D::Digraph)
  result = DigraphWrap.SubdigraphHomeomorphicToK33(GapObj(D))
  @req result != GAP.Globals.fail "subdigraph_homeomorphic_to_k33: the digraph has no subdigraph homeomorphic to K33"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    subdigraph_homeomorphic_to_k4(D::Digraph) -> Vector{Vector{Int}}

Return the immutable list of lists of out-neighbours of a subdigraph of `D`
which is homeomorphic to the complete graph with 4 vertices. Throw an error
if `D` has no such subdigraph.

# Examples
```jldoctest
julia> subdigraph_homeomorphic_to_k4(complete_digraph(5)) isa Vector{Vector{Int}}
true
```
"""
function subdigraph_homeomorphic_to_k4(D::Digraph)
  result = DigraphWrap.SubdigraphHomeomorphicToK4(GapObj(D))
  @req result != GAP.Globals.fail "subdigraph_homeomorphic_to_k4: the digraph has no subdigraph homeomorphic to K4"
  return Vector{Vector{Int}}(result)
end

@doc raw"""
    dual_planar_graph(D::Digraph) -> Digraph

If `D` is a planar digraph, return the dual graph of `D`, i.e. the digraph
which has a vertex for each face of `D` and an edge for each pair of faces
that are separated by an edge of `D`. Throw an error if `D` is not planar.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 4], [1, 4], [2, 3]])
Digraph with 4 vertices, 8 edges

julia> dual_planar_graph(d) isa Digraph
true
```
"""
function dual_planar_graph(D::Digraph)
  result = DigraphWrap.DualPlanarGraph(GapObj(D))
  @req result != GAP.Globals.fail "dual_planar_graph: the digraph is not planar, so its dual does not exist"
  return Digraph(result)
end

@doc raw"""
    digraph_hash(D::Digraph) -> Int

Return a hash value for the digraph `D` based on its structure.

# Examples
```jldoctest
julia> d = digraph([[2], [1]]);

julia> digraph_hash(d)
1541380108904739189
```
"""
function digraph_hash(D::Digraph)
  return Int(DigraphWrap.DigraphHash(GapObj(D))::GapInt)
end

@doc raw"""
    clique_number(D::Digraph) -> Int

Return the clique number of the digraph `D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> clique_number(d)
3
```
"""
function clique_number(D::Digraph)
  return Int(DigraphWrap.CliqueNumber(GapObj(D))::GapInt)
end

@doc raw"""
    has_edge(D::Digraph, s::Int, t::Int) -> Bool

Return `true` if there is at least one edge from `s` to `t` in `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> has_edge(d, 1, 2)
true

julia> has_edge(d, 1, 1)
false
```
"""
function has_edge(D::Digraph, s::Int, t::Int)
  has_vertex(D, s) || return false
  has_vertex(D, t) || return false
  return t in out_neighbours(D)[s]
end

@doc raw"""
    has_vertex(D::Digraph, v::Int) -> Bool

Return `true` if `v` is a valid vertex index in `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> has_vertex(d, 1)
true

julia> has_vertex(d, 3)
false
```
"""
function has_vertex(D::Digraph, v::Int)
  return 1 <= v <= nv(D)
end

@doc raw"""
    graph6_string(D::Digraph) -> String

Return the graph6 representation of the digraph `D`.

# Examples
```jldoctest
julia> d = complete_digraph(2)
Digraph with 2 vertices, 2 edges

julia> graph6_string(d)
"A_"
```
"""
function graph6_string(D::Digraph)
  return String(DigraphWrap.Graph6String(GapObj(D)))
end

@doc raw"""
    sparse6_string(D::Digraph) -> String

Return the sparse6 representation of the digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> sparse6_string(d)
":An"
```
"""
function sparse6_string(D::Digraph)
  return String(DigraphWrap.Sparse6String(GapObj(D)))
end
