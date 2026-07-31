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

Return a topological ordering of the vertices of `D`. The digraph must
be acyclic.

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
  return Vector{Int}(DigraphWrap.DigraphTopologicalSort(GapObj(D)))
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
    digraph_girth(D::Digraph) -> Int

Return the girth (length of the shortest directed cycle) of the digraph `D`.
Returns `-1` if the digraph has no directed cycles.

# Examples
```jldoctest
julia> d = digraph([[2], [3], [1]])
Digraph with 3 vertices, 3 edges

julia> digraph_girth(d)
3
```
"""
function digraph_girth(D::Digraph)
  return Int(DigraphWrap.DigraphGirth(GapObj(D))::GapInt)
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
    digraph_hash(D::Digraph) -> Int

Return a hash value for the digraph `D` based on its structure.

# Examples
```jldoctest
julia> d = digraph([[2], [1]])

julia> digraph_hash(d)
1541380108904739189
```
"""
function digraph_hash(D::Digraph)
  return Int(DigraphWrap.DigraphHash(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_loops(D::Digraph) -> Vector{Tuple{Int,Int}}

Return the list of loops (edges from a vertex to itself) of `D`.

# Examples
```jldoctest
julia> d = digraph([[1, 2], [2]])
Digraph with 2 vertices, 3 edges

julia> digraph_loops(d)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 2)
```
"""
function digraph_loops(D::Digraph)
  return Vector{Tuple{Int,Int}}(DigraphWrap.DigraphLoops(GapObj(D)))
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
    digraph_bicomponents(D::Digraph) -> Union{Nothing, Tuple{Vector{Int}, Vector{Int}}}

Return the bicomponents (bipartition) of the digraph `D` if it is bipartite,
or `nothing` otherwise.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1], [1]])
Digraph with 3 vertices, 4 edges

julia> digraph_bicomponents(d)
2-element Vector{Vector{Int64}}:
 [1]
 [2, 3]
```
"""
function digraph_bicomponents(D::Digraph)
  result = DigraphWrap.DigraphBicomponents(GapObj(D))
  if result == GAP.Globals.fail
    return nothing
  end
  return Tuple{Vector{Int}, Vector{Int}}(result)
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
    chromatic_number(D::Digraph) -> Int

Return the chromatic number of the digraph `D`.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> chromatic_number(d)
3
```
"""
function chromatic_number(D::Digraph)
  return Int(DigraphWrap.ChromaticNumber(GapObj(D))::GapInt)
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
    digraph_shortest_distance(D::Digraph, s::Int, t::Int) -> Int

Return the shortest directed distance from `s` to `t` in `D`.
Returns `-1` if no path exists.

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
  if result == GAP.Globals.fail
    return -1
  end
  return Int(result::GapInt)
end

@doc raw"""
    digraph_path(D::Digraph, s::Int, t::Int) -> Union{Nothing, Tuple{Vector{Int}, Vector{Int}}}

Return a directed path from `s` to `t` in `D` as a list
`[vertices, edges]` where `vertices` lists the vertices along the
path and `edges` contains edge indices (positions in
`OutNeighbours(D)`). Returns GAP `fail` if no path exists.

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
  if result == GAP.Globals.fail
    return nothing
  end
  return Vector{Vector{Int64}}(result)
end

@doc raw"""
    digraph_shortest_path(D::Digraph, s::Int, t::Int) -> Union{Nothing, Tuple{Vector{Int}, Vector{Int}}}

Return a shortest directed path from `s` to `t` in `D` as a list
`[vertices, edges]` (same format as `digraph_path`).
Returns GAP `fail` if no path exists.

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
  if result == GAP.Globals.fail
    return nothing
  end
  return Vector{Vector{Int64}}(result)
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
    digraph_degeneracy(D::Digraph) -> Union{Nothing, Int}

Return the degeneracy of the digraph `D`.

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
  if result == GAP.Globals.fail
    return nothing
  end
  return Int(result::GapInt)
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
    digraph_dijkstra(D::Digraph, s::Int) -> GapObj

Run Dijkstra's algorithm from source vertex `s` on the edge-weighted
digraph `D`.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_dijkstra(d, 1)
2-element Vector{Vector{Int64}}:
 [0, 1, 2]
 [-1, 1, 2]
```
"""
function digraph_dijkstra(D::Digraph, s::Int)
  return Vector{Vector{Int64}}(DigraphWrap.DigraphDijkstra(GapObj(D), s))
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
    digraph_group(D::Digraph) -> PermGroup

Return the group of the Cayley digraph `D` if applicable.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3))

julia> digraph_group(d)
Permutation group of degree 6 and order 6
```
"""
function digraph_group(D::Digraph)
  return PermGroup(DigraphWrap.DigraphGroup(GapObj(D)))
end

@doc raw"""
    group_of_cayley_digraph(D::Digraph) -> PermGroup

Return the group of which `D` is a Cayley digraph, if applicable.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3))

julia> group_of_cayley_digraph(d)
Symmetric group of degree 3
```
"""
function group_of_cayley_digraph(D::Digraph)
  return PermGroup(DigraphWrap.GroupOfCayleyDigraph(GapObj(D)))
end

@doc raw"""
    generators_of_cayley_digraph(D::Digraph) -> GapObj

Return the generating set of the Cayley digraph `D`, if applicable.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3))

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
    digraph_core(D::Digraph) -> Digraph

Return the core of the digraph D.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_core(d)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function digraph_core(D::Digraph)
  return Vector{Int64}(DigraphWrap.DigraphCore(GapObj(D)))
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
    digraph_join_table(D::Digraph) -> Union{Matrix{Int}, Nothing}

Return the join table of the digraph D if D is a join semilattice,
or 
nothing otherwise.

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
  if result == GAP.Globals.fail
    return nothing
  end
  return Matrix{Int}(result)
end

@doc raw"""
    digraph_longest_simple_circuit(D::Digraph) -> Vector{Int}

Return the longest simple circuit in the digraph D.

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
  return Vector{Int}(DigraphWrap.DigraphLongestSimpleCircuit(GapObj(D)))
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
    digraph_meet_table(D::Digraph) -> Union{Matrix{Int}, Nothing}

Return the meet table of the digraph D if D is a meet semilattice,
or 
nothing otherwise.

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
  if result == GAP.Globals.fail
    return nothing
  end
  return Matrix{Int}(result)
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
    digraph_odd_girth(D::Digraph) -> Int

Return the length of the shortest odd directed cycle in D.
Returns -1 if the digraph has no odd directed cycles.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], [1]])
Digraph with 3 vertices, 4 edges

julia> digraph_odd_girth(d)
3
```
"""
function digraph_odd_girth(D::Digraph)
  return Int(DigraphWrap.DigraphOddGirth(GapObj(D))::GapInt)
end

@doc raw"""
    digraph_undirected_girth(D::Digraph) -> Int

Return the girth of the underlying undirected graph of D.
Returns -1 if the underlying undirected graph has no cycles.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> GAP.Globals.DigraphUndirectedGirth(d.X)
3
```
"""
function digraph_undirected_girth(D::Digraph)
  return Int(DigraphWrap.DigraphUndirectedGirth(GapObj(D))::GapInt)
end

@doc raw"""
    edge_weights(D::Digraph) -> Vector{Vector}

Return the edge weights of the edge-weighted digraph D.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2], [1]], [[5], [10]])
Digraph with 2 vertices, 2 edges

julia> edge_weights(d)
2-element Vector{Vector{Int64}}:
 [5]
 [10]
```
"""
function edge_weights(D::Digraph)
  return Vector{Vector}(DigraphWrap.EdgeWeights(GapObj(D)))
end

@doc raw"""
    edge_weighted_digraph_total_weight(D::Digraph) -> GapObj

Return the sum of the edge weights of the edge-weighted digraph D.

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

@doc raw"""
    edge_weighted_digraph_minimum_spanning_tree(D::Digraph) -> Digraph

Return a minimum spanning tree of the edge-weighted digraph D.

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

@doc raw"""
    edge_weighted_digraph_shortest_paths(D::Digraph) -> GapObj

Return the shortest paths data for the edge-weighted digraph D.

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
    hamiltonian_path(D::Digraph) -> Union{Vector{Int}, Nothing}

Return a Hamiltonian path in the digraph D, or 
nothing if none exists.

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
  if result == GAP.Globals.fail
    return nothing
  end
  return Vector{Int}(result)
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
    non_upper_semimodular_pair(D::Digraph) -> Union{Tuple{Int,Int}, Nothing}

Return a pair of vertices witnessing that D is not an upper semimodular
lattice, or 
nothing if no such pair exists.

# Examples
```jldoctest
julia> d = digraph([Vector(1:4), [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> non_upper_semimodular_pair(d)
ERROR: ArgumentError: No such pair exists (meaning that D is an upper semimodular lattice
Stacktrace:
 [1] macro expansion
   @ ~/.julia/packages/AbstractAlgebra/Y7um3/src/Assertions.jl:602 [inlined]
 [2] non_upper_semimodular_pair(D::Digraph)
   @ Oscar /mnt/d/Desktop/Works/Oscar.jl/experimental/Digraphs/src/Attributes.jl:1550
 [3] top-level scope
   @ REPL[81]:1
```
"""
function non_upper_semimodular_pair(D::Digraph)
  result = DigraphWrap.NonUpperSemimodularPair(GapObj(D))
  @req result != GAP.Globals.fail "No such pair exists (meaning that D is an upper semimodular lattice"
  return Vector{Int,Int}(result)
end

@doc raw"""
    non_lower_semimodular_pair(D::Digraph) -> Union{Tuple{Int,Int}, Nothing}

Return a pair of vertices witnessing that D is not a lower semimodular
lattice, or 
nothing if no such pair exists.

# Examples
```jldoctest
julia> d = digraph([Vector(1:4), [2, 4], [3, 4], [4]])
Digraph with 4 vertices, 9 edges

julia> non_lower_semimodular_pair(d)
ERROR: ArgumentError: No such pair exists (meaning that D is an upper semimodular lattice
Stacktrace:
 [1] macro expansion
   @ ~/.julia/packages/AbstractAlgebra/Y7um3/src/Assertions.jl:602 [inlined]
 [2] non_lower_semimodular_pair(D::Digraph)
   @ Oscar /mnt/d/Desktop/Works/Oscar.jl/experimental/Digraphs/src/Attributes.jl:1572
 [3] top-level scope
   @ REPL[82]:1
```
"""
function non_lower_semimodular_pair(D::Digraph)
  result = DigraphWrap.NonLowerSemimodularPair(GapObj(D))
  @req result != GAP.Globals.fail "No such pair exists (meaning that D is an upper semimodular lattice"
  return Vector{Int,Int}(result)
end

@doc raw"""
    strong_orientation(D::Digraph) -> Digraph

Return a strong orientation of the digraph D, if D is bridgeless.

# Examples
```jldoctest
julia> d = cycle_digraph(3)
Digraph with 3 vertices, 3 edges

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
    digraph_maximum_flow(D::Digraph, s::Int, t::Int) -> GapObj

Return the maximum flow from vertex s to vertex 	 in the edge-weighted
digraph D as a GAP record.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2, 2], [3], []], [[3, 2], [1], []])
Digraph with 3 vertices, 3 edges

julia> digraph_maximum_flow(d, 1, 3)
GAP: rec( edges := [ [ 1, 2 ], [ 1, 2 ], [ 2, 3 ] ], flow := [ 2, 1, 1 ], value := 2 )
```
"""
function digraph_maximum_flow(D::Digraph, s::Int, t::Int)
  return DigraphWrap.DigraphMaximumFlow(GapObj(D), s, t)
end

@doc raw"""
    digraph_minimum_cut(D::Digraph, s::Int, t::Int) -> Vector{Vector{Int}}

Return the minimum cut between vertices s and 	 in the edge-weighted
digraph D as a partition of the vertices.

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

@doc raw"""
    digraph_minimum_cut_set(D::Digraph, s::Int, t::Int) -> Vector{Tuple{Int,Int}}

Return the set of edges in a minimum cut between s and 	 in the
edge-weighted digraph D.

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

@doc raw"""
    edge_weighted_digraph_shortest_path(D::Digraph, s::Int, t::Int) -> GapObj

Return the shortest path from s to 	 in the edge-weighted digraph D.

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

@doc raw"""
    edge_weighted_digraph_shortest_paths(D::Digraph, source::Int) -> GapObj

Return the shortest paths from source to every other vertex in the
edge-weighted digraph D.

# Examples
```jldoctest
julia> d = edge_weighted_digraph([[2, 3], [4], [4], []], [[5, 1], [6], [11], []])
Digraph with 4 vertices, 4 edges

julia> edge_weighted_digraph_shortest_paths(d, 1)
GAP: rec( distances := [ 0, 5, 1, 11 ], parents := [ fail, 1, 1, 2 ] )
```
"""
function edge_weighted_digraph_shortest_paths(D::Digraph, source::Int)
  return DigraphWrap.EdgeWeightedDigraphShortestPaths(GapObj(D), source)
end

