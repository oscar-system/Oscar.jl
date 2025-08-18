import Oscar: Polyhedron, Polymake, pm_object
import Oscar.Polymake: Directed, Undirected

function pm_object(G::Graph{T}) where {T <: Union{Directed, Undirected}}
  return G.pm_graph
end

_directed_component(G::MixedGraph) = G.directed_component
@doc raw"""
    directed_component(G::MixedGraph)

Return a copy of the directed graph that is the subgraph of the graph `G`.

# Examples
```jldoctest
julia> mg = graph_from_edges(Mixed, [[1,3],[3,5]],[[4,5],[2,4],[2,3]])
Mixed graph with 5 nodes and the following
Directed edges:
(1, 3)(3, 5)
Undirected edges:
(3, 2)(4, 2)(5, 4)

julia> directed_component(mg)
Directed graph with 5 nodes and the following edges:
(1, 3)(3, 5)
```
"""
directed_component(G::MixedGraph) = deepcopy(_directed_component(G))


_undirected_component(G::MixedGraph) = G.undirected_component
@doc raw"""
    undirected_component(G::MixedGraph)

Return a copy of the directed graph that is the subgraph of the graph `G`.

# Examples
```jldoctest
julia> mg = graph_from_edges(Mixed, [[1,3],[3,5]],[[4,5],[2,4],[2,3]])
Mixed graph with 5 nodes and the following
Directed edges:
(1, 3)(3, 5)
Undirected edges:
(3, 2)(4, 2)(5, 4)

julia> undirected_component(mg)
Undirected graph with 5 nodes and the following edges:
(3, 2)(4, 2)(5, 4)
```
"""
undirected_component(G::MixedGraph) = deepcopy(_undirected_component(G))

################################################################################
################################################################################
##  Constructing and modifying
################################################################################
################################################################################

@doc raw"""
    graph(::Type{T}, nverts::Int64) where {T <: Union{Directed, Mixed, Undirected}}

Construct a graph on `nverts` vertices and no edges. `T` indicates whether the
graph should be `Directed`, `Undirected` or `Mixed`.

# Examples
Make a directed graph with 5 vertices and print the number of nodes and edges.
```jldoctest
julia> g = graph(Directed, 5);

julia> n_vertices(g)
5

julia> n_edges(g)
0
```
"""
function graph(::Type{T}, nverts::Int64) where {T <: Union{Directed, Undirected}}
  return Graph{T}(nverts)
end

function graph(::Type{Mixed}, nverts::Int64)
  return MixedGraph(nverts)
end

@doc raw"""
    graph_from_adjacency_matrix(::Type{T}, G) where {T <:Union{Directed, Undirected}}

Return the graph with adjacency matrix `G`.

This means that the nodes ``i, j`` are connected by an edge
if and only if ``G_{i,j}`` is one.
In the undirected case, it is assumed that ``i > j`` i.e. the upper triangular
part of ``G`` is ignored.

# Examples
```jldoctest
julia> G = ZZ[0 0; 1 0]
[0   0]
[1   0]

julia> graph_from_adjacency_matrix(Directed, G)
Directed graph with 2 nodes and the following edges:
(2, 1)

julia> graph_from_adjacency_matrix(Undirected, G)
Undirected graph with 2 nodes and the following edges:
(2, 1)

```
"""
graph_from_adjacency_matrix(::Type, G::Union{MatElem, Matrix})

function graph_from_adjacency_matrix(::Type{T}, G::Union{MatElem, Matrix}) where {T <: Union{Directed, Undirected}}
  n = nrows(G)
  @req nrows(G)==ncols(G) "not a square matrix"
  g = graph(T, n)
  for i in 1:n
    for j in 1:(T==Undirected ? i-1 : n)
      if isone(G[i,j])
        add_edge!(g, i, j)
      else
        iszero(G[i,j]) || error("not an adjacency matrix")
      end
    end
  end
  return g
end


_has_node(G::Graph, node::Int64) = 0 < node <= n_vertices(G)

@doc raw"""
    add_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}
    add_edge!(mg::MixedGraph, ::Type{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}

Add edge `(s,t)` to the graph `g`.
Return `true` if a new edge `(s,t)` was added, `false` otherwise. For `MixedGraph` the second input determines which component to apply the operation to.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2)
true

julia> add_edge!(g, 1, 2)
false

julia> n_edges(g)
1

julia> mg = graph(Mixed, 3)
Mixed graph with 3 nodes and no edges

julia> add_edge!(mg, Directed, 1, 2)
true

julia> add_edge!(mg, Directed, 1, 2)
false

julia> n_edges(mg)
1

julia> add_edge!(mg, Undirected, 1, 3)
true

julia> n_edges(mg)
2
```
"""
function add_edge!(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
  _has_node(g, source) && _has_node(g, target) || return false
  old_nedges = n_edges(g)
  Polymake._add_edge(pm_object(g), source-1, target-1)
  return n_edges(g) == old_nedges + 1
end

add_edge!(mg::MixedGraph, ::Type{Directed}, source::Int64, target::Int64) = add_edge!(_directed_component(mg), source, target)
add_edge!(mg::MixedGraph, ::Type{Undirected}, source::Int64, target::Int64) = add_edge!(_undirected_component(mg), source, target)

@doc raw"""
    rem_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}
    rem_edge!(g::Graph{T}, e::Edge) where {T <: Union{Directed, Undirected}}
    rem_edge!(g::MixedGraph, ::Type{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}
    rem_edge!(g::MixedGraph, ::Type{T}, e::Edge) where {T <: Union{Directed, Undirected}}

Remove edge `(s,t)` from the graph `g`.
Return `true` if there was an edge from `s` to `t` and it got removed, `false`
otherwise. For `MixedGraph` the second input determines which component to apply the operation to.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2)
true

julia> n_edges(g)
1

julia> rem_edge!(g, 1, 2)
true

julia> n_edges(g)
0

julia> mg = graph(Mixed, 2);

julia> add_edge!(mg, Directed, 1, 2)
true

julia> n_edges(mg)
1

julia> rem_edge!(mg, Directed, 1, 2)
true

julia> add_edge!(mg, Undirected, 1, 2)
true

julia> n_edges(mg)
1

julia> rem_edge!(mg, Undirected, 1, 2)
true

julia> n_edges(mg)
0
```
"""
function rem_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}
  has_edge(g, s, t) || return false
  old_nedges = n_edges(g)
  Polymake._rem_edge(pm_object(g), s-1, t-1)
  return n_edges(g) == old_nedges - 1
end

rem_edge!(mg::MixedGraph, ::Type{Directed}, s::Int64, t::Int64) = rem_edge!(_directed_component(mg), s, t)
rem_edge!(mg::MixedGraph, ::Type{Undirected}, s::Int64, t::Int64) = rem_edge!(_undirected_component(mg), s, t)

@doc raw"""
    add_vertex!(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    add_vertex!(mg::MixedGraph)

Add a vertex to the graph `g`. Return `true` if there a new vertex was actually
added.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> n_vertices(g)
2

julia> add_vertex!(g)
true

julia> n_vertices(g)
3
```
"""
function add_vertex!(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    old_nvertices = n_vertices(g)
    Polymake._add_vertex(pmg)
    return n_vertices(g) - 1 == old_nvertices
end

function add_vertex!(mg::MixedGraph)
  add_vertex!(_directed_component(mg)) || return false
  vertex_added = add_vertex!(_undirected_component(mg))

  if !vertex_added
    @assert rem_vertex!(_directed_component(mg))
    return false
  end
  return true
end

@doc raw"""
    rem_vertex!(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    rem_vertex!(mg::MixedGraph, v::Int64)

Remove the vertex `v` from the graph `g`. Return `true` if node `v` existed and
was actually removed, `false` otherwise.
Please note that this will shift the indices of the vertices with index larger than `v`,
but it will preserve the vertex ordering.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> n_vertices(g)
2

julia> rem_vertex!(g, 1)
true

julia> n_vertices(g)
1
```
"""
function rem_vertex!(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
  _has_node(g, v) || return false
  pmg = pm_object(g)
  old_nvertices = n_vertices(g)
  result = Polymake._rem_vertex(pmg, v-1)
  Polymake._squeeze(pmg)
  return n_vertices(g) + 1 == old_nvertices
end

function rem_vertex!(mg::MixedGraph, v::Int64)
  return rem_vertex!(_directed_component(mg), v) && rem_vertex!(_undirected_component(mg), v)
end

@doc raw"""
    rem_vertices!(g::Graph{T}, a::AbstractArray{Int64}) where {T <: Union{Directed, Undirected}}
    rem_vertices!(mg::MixedGraph, a::AbstractVector)

Remove the vertices in `a` from the graph `g`. Return `true` if at least one vertex was removed.
Please note that this will shift the indices of some of the remaining vertices, but it will preserve the vertex ordering.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> n_vertices(g)
2

julia> rem_vertices!(g, [1, 2])
true

julia> n_vertices(g)
0
```
"""
function rem_vertices!(g::Graph{T}, a::AbstractVector{Int64}) where {T <: Union{Directed, Undirected}}
  pmg = pm_object(g)
  old_nvertices = n_vertices(g)
  for v in a
    0 < v <= old_nvertices && Polymake._rem_vertex(pmg, v-1)
  end
  Polymake._squeeze(pmg)
  return n_vertices(g) < old_nvertices
end

function rem_vertices!(mg::MixedGraph, a::AbstractVector)
  return rem_vertices!(_directed_component(mg), a) && rem_vertices!(_undirected_component(mg), a)
end

@doc raw"""
    add_vertices!(g::Graph{T}, n::Int64) where {T <: Union{Directed, Undirected}}
    add_vertices!(mg::MixedGraph, n::Int64)

Add a `n` new vertices to the graph `g`. Return the number of vertices that
were actually added to the graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> n_vertices(g)
2

julia> add_vertices!(g, 5);

julia> n_vertices(g)
7
```
"""
function add_vertices!(g::Graph{T}, n::Int64) where {T <: Union{Directed, Undirected}}
  return count(_->add_vertex!(g), 1:n)
end

function add_vertices!(mg::MixedGraph, n::Int64)
  return count(_->add_vertex!(mg), 1:n)
end

################################################################################
################################################################################
##  Edges
################################################################################
################################################################################
struct Edge
    source::Int64
    target::Int64
end


@doc raw"""
    src(e::Edge)

Return the source of an edge.

# Examples
```jldoctest
julia> g = complete_graph(2);

julia> E = collect(edges(g));

julia> e = E[1]
Edge(2, 1)

julia> src(e)
2
```
"""
function src(e::Edge)
    return e.source
end


@doc raw"""
    dst(e::Edge)

Return the destination of an edge.

# Examples
```jldoctest
julia> g = complete_graph(2);

julia> E = collect(edges(g));

julia> e = E[1]
Edge(2, 1)

julia> dst(e)
1
```
"""
function dst(e::Edge)
    return e.target
end

Vector{Int}(e::Edge) = [src(e), dst(e)]

Base.isless(a::Edge, b::Edge) = Base.isless(Vector{Int}(a), Vector{Int}(b))

Base.in(i::Int, a::Edge) = (i==src(a) || i==dst(a))

rem_edge!(g::Graph{T}, e::Edge) where {T <: Union{Directed, Undirected}} =
  rem_edge!(g, src(e), dst(e))
rem_edge!(mg::MixedGraph, ::Type{T}, e::Edge) where {T <: Union{Directed, Undirected}} =
  rem_edge!(g, T, src(e), dst(e))

@doc raw"""
    reverse(e::Edge)

Return the edge in the opposite direction of the edge `e`.

# Examples
```jldoctest
julia> g = complete_graph(2);

julia> E = collect(edges(g));

julia> e = E[1]
Edge(2, 1)

julia> reverse(e)
Edge(1, 2)
```
"""
function reverse(e::Edge)
    return Edge(dst(e), src(e))
end

mutable struct EdgeIterator
    pm_itr::Polymake.GraphEdgeIterator{T} where {T <: Union{Directed, Undirected}}
    l::Int64
end
Base.length(eitr::EdgeIterator) = eitr.l
Base.eltype(::Type{EdgeIterator}) = Edge

function Base.iterate(eitr::EdgeIterator, index = 1)
    if eitr.l == 0 || Polymake.isdone(eitr.pm_itr)
        return nothing
    else
        e = Polymake.get_element(eitr.pm_itr)
        s = Polymake.first(e)
        t = Polymake.last(e)
        edge = Edge(s+1, t+1)
        Polymake.increment(eitr.pm_itr)
        eitr.l -= 1
        return (edge, index+1)
    end
end

################################################################################
################################################################################
##  Accessing properties
################################################################################
################################################################################
@doc raw"""
    n_vertices(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the number of vertices of a graph.

# Examples
The edge graph of the cube has eight vertices, just like the cube itself.
```jldoctest
julia> c = cube(3);

julia> g = vertex_edge_graph(c);

julia> n_vertices(g)
8
```
"""
function n_vertices(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return Polymake.nv(pm_object(g))
end

function n_vertices(g::MixedGraph)
  return n_vertices(_directed_component(g))
end


@doc raw"""
    vertices(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    vertices(g::MixedGraph)

Return the vertex indices of a graph.

# Examples
The edge graph of the cube has eight vertices, numbered 1 to 8.
```jldoctest
julia> c = cube(3);

julia> g = vertex_edge_graph(c);

julia> vertices(g)
1:8
```
"""
function vertices(g::Graph{T}) where {T <: Union{Directed, Undirected}}
  return 1:n_vertices(g)
end

function vertices(g::MixedGraph)
  return 1:n_vertices(g)
end

@doc raw"""
    n_edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    n_edges(g::MixedGraph)

Return the number of edges of a graph.

# Examples
The edge graph of the cube has 12 edges just like the cube itself.
```jldoctest
julia> c = cube(3);

julia> g = vertex_edge_graph(c);

julia> n_edges(g)
12
```
"""
function n_edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return Polymake.ne(pm_object(g))
end

function n_edges(g::MixedGraph)
  return sum(Polymake.ne.(pm_object.([_directed_component(g),
                                      _undirected_component(g)])))
end

@doc raw"""
    edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    edges(mg::MixedGraph)

Return an iterator over the edges of the graph `g`.

# Examples
A triangle has three edges.
```jldoctest
julia> triangle = simplex(2);

julia> g = vertex_edge_graph(triangle);

julia> collect(edges(g))
3-element Vector{Edge}:
 Edge(2, 1)
 Edge(3, 1)
 Edge(3, 2)

julia> mg = graph_from_edges(Mixed, [[1,3],[3,5]],[[4,5],[2,4],[2,3]])
Mixed graph with 5 nodes and the following
Directed edges:
(1, 3)(3, 5)
Undirected edges:
(3, 2)(4, 2)(5, 4)

julia> collect(edges(mg, Directed))
2-element Vector{Edge}:
 Edge(1, 3)
 Edge(3, 5)

```
"""
function edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return EdgeIterator(Polymake.edgeiterator(pm_object(g)), n_edges(g))
end

edges(g::MixedGraph, ::Type{Directed}) = edges(_directed_component(g))
edges(g::MixedGraph, ::Type{Undirected}) = edges(_undirected_component(g))

@doc raw"""
    has_edge(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
    has_edge(g::MixedGraph, ::Type{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}

Check for an edge in a graph.

# Examples
Check for the edge $1\to 2$ in the edge graph of a triangle.
```jldoctest
julia> triangle = simplex(2);

julia> g = vertex_edge_graph(triangle);

julia> has_edge(g, 1, 2)
true

julia> mg = graph_from_edges(Mixed, [[1,3],[3,5]],[[4,5],[2,4],[2,3]])
Mixed graph with 5 nodes and the following
Directed edges:
(1, 3)(3, 5)
Undirected edges:
(3, 2)(4, 2)(5, 4)

julia> has_edge(mg, Directed, 1, 3)
true

julia> has_edge(mg, Undirected, 1, 3)
false
```
"""
function has_edge(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    return Polymake._has_edge(pmg, source-1, target-1)
end
function has_edge(g::Graph{T}, e::Edge) where {T <: Union{Directed, Undirected}}
    return has_edge(g, src(e), dst(e))
end

has_edge(mg::MixedGraph, ::Type{Directed}, source::Int64, target::Int64) = has_edge(_directed_component(mg), source, target)
has_edge(mg::MixedGraph, ::Type{Undirected}, source::Int64, target::Int64) = has_edge(_undirected_component(mg), source, target)

function has_edge(g::MixedGraph, e::Edge)
    return has_edge(g, T, src(e), dst(e))
end

@doc raw"""
    has_vertex(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    has_vertex(g::MixedGraph, v::Int64)
Check for a vertex in a graph.

# Examples
The edge graph of a triangle only has 3 vertices.
```jldoctest
julia> triangle = simplex(2);

julia> g = vertex_edge_graph(triangle);

julia> has_vertex(g, 1)
true

julia> has_vertex(g, 4)
false
```
"""
function has_vertex(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
  pmg = pm_object(g)
  return Polymake._has_vertex(pmg, v-1)
end

#each component has same vertices so it is enough to check one
has_vertex(g::MixedGraph, v::Int) = has_vertex(_directed_component(g), v) 

@doc raw"""
    neighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}

Return the neighboring vertices of a vertex `v` in a graph `g`. If the graph is
directed, the neighbors reachable via outgoing edges are returned.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> neighbors(g, 3)
1-element Vector{Int64}:
 4
```
"""
function neighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g);
    result = Polymake._outneighbors(pmg, v-1)
    return [x+1 for x in result]
end

@doc raw"""
    degree(g::Graph{T} [, v::Int64]) where {T <: Union{Directed, Undirected}}

Return the degree of the vertex `v` in the graph `g`. If `v` is
missing, return the list of degrees of all vertices. If the graph is
directed, only neighbors reachable via outgoing edges are counted.

See also [`indegree`](@ref) and [`outdegree`](@ref) for directed graphs.

# Examples
```jldoctest
julia> g = vertex_edge_graph(icosahedron());

julia> degree(g, 1)
5
```
"""
degree(g::Graph, v::Int64) = length(neighbors(g, v))

@doc raw"""
    indegree(g::Graph{Directed} [, v::Int64])

Return the indegree of the vertex `v` in the directed graph `g`. If `v` is
missing, return the list of indegrees of all vertices.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> indegree(g, 1)
0

julia> indegree(g)
5-element Vector{Int64}:
 0
 0
 1
 1
 0
```
"""
indegree(g::Graph{Directed}, v::Int64) = length(inneighbors(g, v))

@doc raw"""
    outdegree(g::Graph{Directed} [, v::Int64])

Return the outdegree of the vertex `v` in the directed graph `g`. If `v` is
missing, return the list of outdegrees of all vertices.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> outdegree(g, 1)
1

julia> outdegree(g)
5-element Vector{Int64}:
 1
 0
 1
 0
 0
```
"""
outdegree(g::Graph{Directed}, v::Int64) = length(outneighbors(g, v))

degree(g::Graph) = [ length(neighbors(g, v)) for v in 1:n_vertices(g) ]
indegree(g::Graph{Directed}) = [ length(inneighbors(g, v)) for v in 1:n_vertices(g) ]
outdegree(g::Graph{Directed}) = [ length(outneighbors(g, v)) for v in 1:n_vertices(g) ]


@doc raw"""
    inneighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}

Return the vertices of a graph `g` that have an edge going towards `v`. For an
undirected graph, all neighboring vertices are returned.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> inneighbors(g, 3)
1-element Vector{Int64}:
 1

julia> inneighbors(g, 1)
Int64[]
```
"""
function inneighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g);
    result = Polymake._inneighbors(pmg, v-1)
    return [x+1 for x in result]
end


@doc raw"""
    outneighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}

Return the vertices of a graph `g` that are target of an edge coming from `v`.
For an undirected graph, all neighboring vertices are returned.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> outneighbors(g, 3)
1-element Vector{Int64}:
 4

julia> outneighbors(g, 4)
Int64[]
```
"""
function outneighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g);
    result = Polymake._outneighbors(pmg, v-1)
    return [x+1 for x in result]
end


@doc raw"""
    all_neighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}

Return all vertices of a graph `g` that are connected to the vertex `v` via an
edge, independent of the edge direction.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> all_neighbors(g, 3)
2-element Vector{Int64}:
 1
 4

julia> all_neighbors(g, 4)
1-element Vector{Int64}:
 3
```
"""
function all_neighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g);
    result = union(Polymake._inneighbors(pmg, v-1), Polymake._outneighbors(pmg, v-1))
    return [x+1 for x in result]
end

@doc raw"""
    incidence_matrix(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return an unsigned (boolean) incidence matrix representing a graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g, 1, 3);

julia> add_edge!(g, 3, 4);

julia> incidence_matrix(g)
5×2 IncidenceMatrix
 [1]
 []
 [1, 2]
 [2]
 []
```
"""
function incidence_matrix(g::Graph{T}) where {T <: Union{Directed, Undirected}}
  IncidenceMatrix(Polymake.graph.incidence_matrix(pm_object(g)))
end

@doc raw"""
    signed_incidence_matrix(g::Graph)

Return a signed incidence matrix representing a graph `g`.  If `g` is directed, sources will have sign `-1` and targest will have sign `+1`.  If `g` is undirected, vertices of larger index will have sign `-1` and vertices of smaller index will have sign `+1`.

# Examples
```jldoctest
julia> g = Graph{Directed}(5);

julia> add_edge!(g,1,2); add_edge!(g,2,3); add_edge!(g,3,4); add_edge!(g,4,5); add_edge!(g,5,1);

julia> signed_incidence_matrix(g)
5×5 Matrix{Int64}:
 -1   0   0   0   1
  1  -1   0   0   0
  0   1  -1   0   0
  0   0   1  -1   0
  0   0   0   1  -1

julia> g = Graph{Undirected}(5);

julia> add_edge!(g,1,2); add_edge!(g,2,3); add_edge!(g,3,4); add_edge!(g,4,5); add_edge!(g,5,1);

julia> signed_incidence_matrix(g)
5×5 Matrix{Int64}:
  1   0   0   1   0
 -1   1   0   0   0
  0  -1   1   0   0
  0   0  -1   0   1
  0   0   0  -1  -1

```
"""
signed_incidence_matrix(g::Graph) = convert(Matrix{Int}, Polymake.graph.signed_incidence_matrix(pm_object(g)))

# EdgeMap getters and setters for Edges
function Base.getindex(M::Polymake.EdgeMap{TK, TV}, e::Edge) where {TK, TV}
  return M[src(e), dst(e)]
end

function Base.setindex!(M::Polymake.EdgeMap{T, TV}, val, e::Edge) where {T, TV}
  M[src(e), dst(e)] = val
  return val
end

_has_vertex_map(GM::GraphMap) = !isnothing(GM.vertex_map)
_has_edge_map(GM::GraphMap) = !isnothing(GM.edge_map)

function Base.getindex(GM::GraphMap, i::Int)
  @req _has_vertex_map(GM) "Graph map not defined for vertices"
  @req _has_node(GM.graph, i) "Graph doesn't have vertex $i"
  return GM.vertex_map[i]
end

function Base.getindex(GM::GraphMap, i::Int, j::Int)
  @req has_edge(GM.graph, i, j) "Graph doesn't have edge ($i, $j) "
  @req _has_edge_map(GM) "Graph map not defined on edges"
  return GM.edge_map[i, j]
end

function Base.getindex(GM::GraphMap, e::Edge)
  return GM[src(e), dst(e)]
end

function Base.setindex!(GM::GraphMap, val, i::Int)
  @req _has_node(GM.graph, i) "Graph doesn't have vertex $i"
  GM.vertex_map[i] = val
  return val
end

function Base.setindex!(GM::GraphMap, val, i::Int, j::Int)
  @req has_edge(GM.graph, i, j) "Graph doesn't have edge ($i, $j) "
  GM.edge_map[i, j] = val
  return val
end

function Base.setindex!(GM::GraphMap, val, indices::Tuple{Int, Int})
  GM[indices[1], indices[2]] = val
  return val
end

function Base.setindex!(GM::GraphMap, val, e::Edge)
  GM[src(e), dst(e)] = val
  return val
end

function Base.getproperty(G::Graph, p::Symbol)
  hasfield(Graph, p) && return getfield(G, p)
  @req has_attribute(G, p) "$G doesn't have a labeling $p"
  return get_attribute(G, p)
end

@doc raw"""
    labelings(G::Graph)

Return the names of all labelings of a graph `G` as a `Vector{Symbol}`.

# Examples
```jldoctest
julia> G = graph_from_labeled_edges(Directed, Dict((1, 2) => 1, (2, 3) => 4); name=:color)
Directed graph with 3 nodes and the following labeling(s):
label: color
(1, 2) -> 1
(2, 3) -> 4

julia> labelings(G)
1-element Vector{Symbol}:
 :color

```
"""
function labelings(G::AbstractGraph)
  attrs = AbstractAlgebra._get_attributes(G)
  isnothing(attrs) && return Symbol[]
  return [k for (k, v ) in attrs if v isa GraphMap]
end

function _graph_maps(G::Graph)
  labels = tuple(labelings(G)...)
  isempty(labels) && return Dict{Symbol, GraphMap}()
  return Dict(l => getproperty(G, l) for l in labels)
end

################################################################################
################################################################################
##  Higher order algorithms
################################################################################
################################################################################
@doc raw"""
    automorphism_group_generators(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return generators of the automorphism group of the graph `g`.

# Examples
```jldoctest
julia> g = complete_graph(4);

julia> automorphism_group_generators(g)
3-element Vector{PermGroupElem}:
 (3,4)
 (2,3)
 (1,2)
```
"""
function automorphism_group_generators(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g);
    result = Polymake.graph.automorphisms(pmg)
    return _pm_arr_arr_to_group_generators(result, n_vertices(g))
end


@doc raw"""
    automorphism_group(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the automorphism group of the graph `g`.

# Examples
```jldoctest
julia> g = complete_graph(4);

julia> automorphism_group(g)
Permutation group of degree 4
```
"""
function automorphism_group(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return _gens_to_group(automorphism_group_generators(g))
end


@doc raw"""
    shortest_path_dijkstra(g::Graph{T}, s::Int64, t::Int64; reverse::Bool=false) where {T <: Union{Directed, Undirected}}

Compute the shortest path between two vertices in a graph using Dijkstra's
algorithm. All edges are set to have a length of 1. The optional parameter
indicates whether the edges should be considered reversed.

# Examples
```jldoctest
julia> g = Graph{Directed}(3);

julia> add_edge!(g, 1, 2);

julia> add_edge!(g, 2, 3);

julia> add_edge!(g, 3, 1);

julia> shortest_path_dijkstra(g, 3, 1)
2-element Vector{Int64}:
 3
 1

julia> shortest_path_dijkstra(g, 1, 3)
3-element Vector{Int64}:
 1
 2
 3

julia> shortest_path_dijkstra(g, 3, 1; reverse=true)
3-element Vector{Int64}:
 3
 2
 1
```
"""
function shortest_path_dijkstra(g::Graph{T}, s::Int64, t::Int64; reverse::Bool=false) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    em = Polymake.EdgeMap{T, Int64}(pmg)
    for e in edges(g)
        Polymake._set_entry(em, src(e)-1, dst(e)-1, 1)
    end
    result = Polymake._shortest_path_dijkstra(pmg, em, s-1, t-1, !reverse)
    return Polymake.to_one_based_indexing(result)
end

@doc raw"""
    connectivity(g::Graph{Undirected})

Return the connectivity of the undirected graph `g`.

# Examples
```jldoctest
julia> g = complete_graph(3);

julia> connectivity(g)
2

julia> rem_edge!(g, 2, 3);

julia> connectivity(g)
1

julia> rem_edge!(g, 1, 3);

julia> connectivity(g)
0
```
"""
connectivity(g::Graph{Undirected}) = Polymake.graph.connectivity(g)::Int

@doc raw"""
    is_connected(g::Graph{Undirected})

Check if the undirected graph `g` is connected.

# Examples
```jldoctest
julia> g = Graph{Undirected}(3);

julia> is_connected(g)
false

julia> add_edge!(g, 1, 2);

julia> add_edge!(g, 2, 3);

julia> is_connected(g)
true
```
"""
is_connected(g::Graph{Undirected}) = Polymake.call_function(:graph, :is_connected, pm_object(g))::Bool

@doc raw"""
    connected_components(g::Graph{Undirected})

Return the connected components of an undirected graph `g`.

# Examples
```jldoctest
julia> g = Graph{Undirected}(2);

julia> add_edge!(g, 1, 2);

julia> connected_components(g)
1-element Vector{Vector{Int64}}:
 [1, 2]
```
"""
function connected_components(g::Graph{Undirected})
    im = Polymake.call_function(:graph, :connected_components, pm_object(g))::IncidenceMatrix
    return [Vector(Polymake.row(im,i)) for i in 1:Polymake.nrows(im)]
end

@doc raw"""
    is_strongly_connected(g::Graph{Directed})

Check if the directed graph `g` is strongly connected.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> is_strongly_connected(g)
false

julia> add_edge!(g, 2, 1);

julia> is_strongly_connected(g)
true
```
"""
is_strongly_connected(g::Graph{Directed}) = Polymake.call_function(:graph, :is_strongly_connected, pm_object(g))::Bool

@doc raw"""
    strongly_connected_components(g::Graph{Directed})

Return the strongly connected components of a directed graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> length(strongly_connected_components(g))
2

julia> add_edge!(g, 2, 1);

julia> strongly_connected_components(g)
1-element Vector{Vector{Int64}}:
 [1, 2]
```
"""
function strongly_connected_components(g::Graph{Directed})
    im = Polymake.call_function(:graph, :strong_components, pm_object(g))::IncidenceMatrix
    return [Vector(Polymake.row(im,i)) for i in 1:Polymake.nrows(im)]
end

@doc raw"""
    is_weakly_connected(g::Graph{Directed})

Check if the directed graph `g` is weakly connected.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> is_weakly_connected(g)
true
```
"""
is_weakly_connected(g::Graph{Directed}) = Polymake.call_function(:graph, :is_weakly_connected, pm_object(g))::Bool

@doc raw"""
    weakly_connected_components(g::Graph{Directed})

Return the weakly connected components of a directed graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> weakly_connected_components(g)
1-element Vector{Vector{Int64}}:
 [1, 2]
```
"""
function weakly_connected_components(g::Graph{Directed})
    im = Polymake.call_function(:graph, :weakly_connected_components, pm_object(g))::IncidenceMatrix
    return [Vector(Polymake.row(im,i)) for i in 1:Polymake.nrows(im)]
end

@doc raw"""
    diameter(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the diameter of a (strongly) connected (di-)graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> weakly_connected_components(g)
1-element Vector{Vector{Int64}}:
 [1, 2]
```
"""
function diameter(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    if T == Directed && !is_strongly_connected(g) ||
       T == Undirected && !is_connected(g)
        throw(ArgumentError("The (di-)graph must be (strongly) connected!"))
    end
    return Polymake.call_function(:graph, :diameter, pm_object(g))::Int
end

@doc raw"""
    is_isomorphic(g1::Graph{T}, g2::Graph{T}) where {T <: Union{Directed, Undirected}}

Check if the graph `g1` is isomorphic to the graph `g2`.

# Examples
```jldoctest
julia> is_isomorphic(vertex_edge_graph(simplex(3)), dual_graph(simplex(3)))
true

julia> is_isomorphic(vertex_edge_graph(cube(3)), dual_graph(cube(3)))
false
```
"""
is_isomorphic(g1::Graph{T}, g2::Graph{T}) where {T <: Union{Directed, Undirected}} = Polymake.graph.isomorphic(pm_object(g1), pm_object(g2))::Bool

@doc raw"""
    is_isomorphic_with_permutation(G1::Graph, G2::Graph) -> Bool, Vector{Int}

Return whether `G1` is isomorphic to `G2` as well as a permutation
of the nodes of `G1` such that both graphs agree.

# Examples
```jldoctest
julia> is_isomorphic_with_permutation(vertex_edge_graph(simplex(3)), dual_graph(simplex(3)))
(true, [1, 2, 3, 4])

```
"""
function is_isomorphic_with_permutation(G1::Graph, G2::Graph)
  f12 = Polymake.graph.find_node_permutation(G1.pm_graph, G2.pm_graph)
  if isnothing(f12)
    return false, Vector{Int}()
  end
  return true, Polymake.to_one_based_indexing(f12)
end

@doc raw"""
    _is_equal_up_to_permutation_with_permutation(A1::MatElem, A2::MatElem) -> Bool, Vector{Int}

Return a permutation `I` such that `A1[I,I] == A2` and whether it exists.

The method assumes that both matrices are symmetric, their diagonal entries
are all equal (and so irrelevant) and the off-diagonal entries are either ``0``
or ``1``. It is assumed that `A1` and `A2` are symmetric and
their upper triangular part is ignored.
"""
function _is_equal_up_to_permutation_with_permutation(A1::MatElem, A2::MatElem)
  g1 = graph_from_adjacency_matrix(Undirected, A1)
  g2 = graph_from_adjacency_matrix(Undirected, A2)
  b, T = is_isomorphic_with_permutation(g1, g2)
  if b
    @assert A1[T, T] == A2
  end
  return b, T
end


################################################################################
################################################################################
##  Standard constructions
################################################################################
################################################################################
@doc raw"""
    vertex_edge_graph(p::Polyhedron)

Return the edge graph of a `Polyhedron`, vertices of the graph correspond to
vertices of the polyhedron, there is an edge between two vertices if the
polyhedron has an edge between the corresponding vertices. The resulting graph
is `Undirected`.
If the polyhedron has lineality, then it has no vertices or bounded edges, so the `vertex_edge_graph` will be the empty graph.
In this case, the keyword argument can be used to consider the polyhedron modulo its lineality space.

# Examples
Construct the edge graph of the cube. Like the cube it has 8 vertices and 12
edges.
```jldoctest
julia> c = cube(3);

julia> g = vertex_edge_graph(c);

julia> n_vertices(g)
8

julia> n_edges(g)
12
```
"""
function vertex_edge_graph(p::Polyhedron; modulo_lineality=false)
  lineality_dim(p) != 0 && !modulo_lineality && return Graph{Undirected}(0)
  og = Graph{Undirected}(pm_object(p).GRAPH.ADJACENCY)
  is_bounded(p) || rem_vertices!(og, _ray_indices(pm_object(p)))
  return og
end

@doc raw"""
    dual_graph(p::Polyhedron)

Return the dual graph of a `Polyhedron`, vertices of the graph correspond to
facets of the polyhedron and there is an edge between two vertices if the
corresponding facets are neighboring, meaning their intersection is a
codimension 2 face of the polyhedron.

For bounded polyhedra containing 0 in the interior this is the same as the edge
graph the polar dual polyhedron.

# Examples
Construct the dual graph of the cube. This is the same as the edge graph of the
octahedron, so it has 6 vertices and 12 edges.
```jldoctest
julia> c = cube(3);

julia> g = dual_graph(c);

julia> n_vertices(g)
6

julia> n_edges(g)
12
```
"""
function dual_graph(p::Polyhedron)
  pop = pm_object(p)
  og = Graph{Undirected}(pop.DUAL_GRAPH.ADJACENCY)
  fai = _facet_at_infinity(pop)
  if fai <= n_vertices(og)
    rem_vertex!(og, fai)
  elseif !is_bounded(p)
    f = facets(IncidenceMatrix, p)
    farface = _ray_indices(pop)
    for e in edges(og)
      if is_subset(intersect(row(f, src(e)), row(f, dst(e))), farface)
        rem_edge!(og, e)
      end
    end
  end
  return og
end

@doc raw"""
    dual_graph(SOP::SubdivisionOfPoints)

Return the dual graph of a `SubdivisionOfPoints`, nodes correspond to the maximal cells,
and there is an edge if two maximal cells share a common face of codimension one.

# Examples
Construct the dual graph of a triangulation of the square; it has a single edge.
```jldoctest
julia> S = subdivision_of_points(vertices(cube(2)), [0,0,0,1])
Subdivision of points in ambient dimension 2

julia> dual_graph(S)
Undirected graph with 2 nodes and the following edges:
(2, 1)
```
"""
dual_graph(SOP::SubdivisionOfPoints) = Graph{Undirected}(pm_object(SOP).POLYHEDRAL_COMPLEX.DUAL_GRAPH.ADJACENCY)


@doc raw"""
    complete_graph(n::Int64)

Assemble the undirected complete graph on `n` nodes.

# Examples
```jldoctest
julia> g = complete_graph(3);

julia> collect(edges(g))
3-element Vector{Edge}:
 Edge(2, 1)
 Edge(3, 1)
 Edge(3, 2)
```
"""
function complete_graph(n::Int64)
    bigobj = Polymake.graph.complete(n)
    return Graph{Undirected}(bigobj.ADJACENCY)
end


@doc raw"""
    complete_bipartite_graph(n::Int64, m::Int64)

Assemble the undirected complete bipartite graph between `n` and `m` nodes.

# Examples
```jldoctest
julia> g = complete_bipartite_graph(2,2);

julia> collect(edges(g))
4-element Vector{Edge}:
 Edge(3, 1)
 Edge(3, 2)
 Edge(4, 1)
 Edge(4, 2)
```
"""
function complete_bipartite_graph(n::Int64, m::Int64)
    bigobj = Polymake.graph.complete_bipartite(n, m)
    return Graph{Undirected}(bigobj.ADJACENCY)
end



@doc raw"""
    visualize(G::Graph{<:Union{Polymake.Directed, Polymake.Undirected}}; backend::Symbol=:threejs, filename::Union{Nothing, String}=nothing, kwargs...)

Visualize a graph, see [`visualize`](@ref Oscar.visualize(::Union{SimplicialComplex, Cone{<:Union{Float64, FieldElem}}, Graph, PolyhedralComplex{<:Union{Float64, FieldElem}}, PolyhedralFan{<:Union{Float64, FieldElem}}, Polyhedron, SubdivisionOfPoints{<:Union{Float64, FieldElem}}})) for details on the keyword arguments.

The `backend` keyword argument allows the user to pick between a Three.js visualization by default, or passing `:tikz` for a TikZ visualization.
The `filename` keyword argument will write visualization code to the `filename` location, this will be html for `:threejs` backend or TikZ code for `:tikz`.
If the graph `G` has a labeling `:color` (see [`label!`](@ref)) then the visualization will use these colors to color the graph.

Possible color labelings include RGB values of the form `"255 0 255"` or `"#ff00ff"`, as well as the following named colors as strings: `polymakeorange`, `polymakegreen`,
`white`, `purple`, `cyan`, `darkolivegreen`, `indianred`, `plum1`, `red`, `lightslategrey`, `yellow`, `orange`, `salmon1`, `azure`, `green`, `gray`, `midnightblue`, `pink`, `magenta`, `blue`, `lavenderblush`, `chocolate1`, `lightgreen`, `black`.

"""
function visualize(G::Graph{T};
                   backend::Symbol=:threejs, filename::Union{Nothing, String}=nothing,
                   kwargs...) where {T <: Union{Directed, Undirected}}
  BG = Polymake.graph.Graph{T}(ADJACENCY=pm_object(G))

  defaults = (;VertexLabels = collect(1:n_vertices(G)))
  if has_attribute(G, :color)
    defaults = merge(defaults,
                     NamedTuple(k => v for (k, v) in
                                  [(:EdgeColor, G.color.edge_map), (:VertexColor, G.color.vertex_map)] if !isnothing(v)))
  end
  
  vBG = Polymake.visual(Polymake.Visual, BG; merge(defaults, kwargs)...)
  if !isnothing(filename)
    Polymake.call_function(Nothing, :graph, backend, vBG; File=filename)
  else
    Polymake.call_function(Nothing, :graph, backend, vBG;)
  end
end


# Some standard polytopes from graphs
@doc raw"""
    fractional_cut_polytope(G::Graph{Undirected})

Construct the fractional cut polytope of the graph $G$.


# Examples
```jldoctest
julia> G = complete_graph(4);

julia> fractional_cut_polytope(G)
Polytope in ambient dimension 6
```
"""
fractional_cut_polytope(G::Graph{Undirected}) = polyhedron(Polymake.polytope.fractional_cut_polytope(pm_object(G)))


@doc raw"""
    fractional_matching_polytope(G::Graph{Undirected})

Construct the fractional matching polytope of the graph $G$.


# Examples
```jldoctest
julia> G = complete_graph(4);

julia> fractional_matching_polytope(G)
Polytope in ambient dimension 6
```
"""
fractional_matching_polytope(G::Graph{Undirected}) = polyhedron(Polymake.polytope.fractional_matching_polytope(pm_object(G)))


################################################################################
################################################################################
##  Printing
################################################################################
################################################################################
_to_string(::Type{Polymake.Directed}) = "Directed"
_to_string(::Type{Polymake.Undirected}) = "Undirected"
_to_string(::Type{Mixed}) = "Mixed"

function Base.show(io::IO, m::MIME"text/plain", G::AbstractGraph{T}) where {T <: GraphTypes}
  if n_edges(G) > 0
    labels = labelings(G)
    if !isempty(labels)
      print(io, "$(_to_string(T)) graph with $(n_vertices(G)) nodes and the following labeling(s):")
      for label in labels
        println(io, "")
        print(io, "label: $label")
        if _has_edge_map(getproperty(G, label))
          for e in edges(G)
            println(io, "")
            print(io, "($(src(e)), $(dst(e))) -> $(getproperty(G, label)[e])")
          end
        end
        if _has_vertex_map(getproperty(G, label))
          for v in 1:n_vertices(G)
            println(io, "")
            print(io, "$v -> $(getproperty(G, label)[v])")
          end
        end
      end
    else
      if T == Mixed
        println(io, "$(_to_string(T)) graph with $(n_vertices(G)) nodes and the following")  # at least one new line is needed
        println(io, "Directed edges:")
        for e in edges(G, Directed)
          print(io, "($(src(e)), $(dst(e)))")
        end
        println(io, "")
        println(io, "Undirected edges:")
        for e in edges(G, Undirected)
          print(io, "($(src(e)), $(dst(e)))")
        end
      else
        println(io, "$(_to_string(T)) graph with $(n_vertices(G)) nodes and the following edges:")  # at least one new line is needed
        for e in edges(G)
          print(io, "($(src(e)), $(dst(e)))")
        end
      end
    end
  else
    print(io, "$(_to_string(T)) graph with $(n_vertices(G)) nodes and no edges")
  end
end

function Base.show(io::IO, G::AbstractGraph{T})  where {T <: GraphTypes}
  if is_terse(io)
    !isempty(labelings(G)) && print(io, "Labeled ")
    print(io, "$(_to_string(T)) graph")
  else
    print(io, "$(_to_string(T)) graph with $(n_vertices(G)) nodes and $(n_edges(G)) edges")
    !isempty(labelings(G)) && print(io, " with labeling(s) $(labelings(G))")
  end
end

function graph_from_edges(::Type{T},
                          edges::Vector{Edge},
                          n_vertices::Int=-1) where {T <: Union{Directed, Undirected}}
  n_needed = maximum(reduce(append!,[[src(e),dst(e)] for e in edges]; init=[0]))
  @req (n_vertices >= n_needed || n_vertices < 0)  "n_vertices must be at least the maximum vertex in the edges"

  g = graph(T, max(n_needed, n_vertices))
  for e in edges
    add_edge!(g, src(e), dst(e))
  end

  return g
end

function graph_from_edges(::Type{T},
                          edges::EdgeIterator,
                          n_vertices::Int=-1) where {T <: Union{Directed, Undirected}}
  return graph_from_edges(T, collect(edges), n_vertices)
end

@doc raw"""
    graph_from_edges(edges::Vector{Vector{Int}})
    graph_from_edges(::Type{T}, edges::Vector{Vector{Int}}, n_vertices::Int=-1) where {T <:Union{Directed, Undirected}}
    graph_from_edges(::Type{Mixed}, directed_edges::Vector{Vector{Int}}, undirected_edges::Vector{Vector{Int}}; n_vertices=-1)
    graph_from_edges(::Type{Mixed}, directed_edges::Vector{Edge}, undirected_edges::Vector{Edge}; n_vertices=-1)


Create a graph from a vector of edges. There is an optional input for number of vertices, `graph_from_edges`  will
ignore any negative integers and throw an error when the input is less than the maximum vertex index in edges.

# Examples
```jldoctest
julia> G = graph_from_edges([[1,3],[3,5],[4,5],[2,4],[2,3]])
Undirected graph with 5 nodes and the following edges:
(3, 1)(3, 2)(4, 2)(5, 3)(5, 4)

julia> G = graph_from_edges(Directed, [[1,3]], 4)
Directed graph with 4 nodes and the following edges:
(1, 3)

julia> g = graph_from_edges(Mixed, [[1,3],[3,5]],[[4,5],[2,4],[2,3]])
Mixed graph with 5 nodes and the following
Directed edges:
(1, 3)(3, 5)
Undirected edges:
(3, 2)(4, 2)(5, 4)

julia> g1 = graph_from_edges(Directed, [[1, 2], [3, 4]])
Directed graph with 4 nodes and the following edges:
(1, 2)(3, 4)

julia> g2 = graph_from_edges(Undirected, [[1, 4], [3, 2]])
Undirected graph with 4 nodes and the following edges:
(3, 2)(4, 1)

julia> graph_from_edges(Mixed, edges(g1), edges(g2))
Mixed graph with 4 nodes and the following
Directed edges:
(1, 2)(3, 4)
Undirected edges:
(3, 2)(4, 1)
```
"""
function graph_from_edges(::Type{T},
                          edges::Vector{S},
                          n_vertices::Int=-1) where {T <: Union{Directed, Undirected}, S <: Union{Vector{Int}, NTuple{2, Int}}}
  return graph_from_edges(T, [Edge(e[1], e[2]) for e in edges], n_vertices)
end

function graph_from_edges(edges::Vector{T},
                          n_vertices::Int=-1) where T <: Union{Vector{Int}, NTuple{2, Int}}
  return graph_from_edges(Undirected, [Edge(e[1], e[2]) for e in edges], n_vertices)
end


function graph_from_edges(::Type{Mixed},
                          directed_edges::Vector{Edge},
                          undirected_edges::Vector{Edge},
                          n_vertices::Int=-1)
  n_needed = maximum(vcat(reduce(append!,[[src(e),dst(e)] for e in directed_edges]; init=[0]),
                          reduce(append!,[[src(e),dst(e)] for e in undirected_edges]; init=[0])))
  @req (n_vertices >= n_needed || n_vertices < 0)  "n_vertices must be at least the maximum vertex in the edges"
  
  g = graph(Mixed, max(n_needed, n_vertices))
  for e in directed_edges
    add_edge!(g, Directed, src(e), dst(e))
  end

  for e in undirected_edges
    add_edge!(g, Undirected, src(e), dst(e))
  end

  return g
end

function graph_from_edges(::Type{Mixed},
                          directed_edges::Vector{S},
                          undirected_edges::Vector{T},
                          n_vertices::Int=-1) where {S, T <: Union{Vector{Int}, NTuple{2, Int}}}
  return graph_from_edges(Mixed,
                          [Edge(e[1], e[2]) for e in directed_edges],
                          [Edge(e[1], e[2]) for e in undirected_edges],
                          n_vertices)
end

function graph_from_edges(::Type{Mixed},
                          directed_edges::EdgeIterator,
                          undirected_edges::EdgeIterator,
                          n_vertices::Int=-1) 
  return graph_from_edges(Mixed, collect(directed_edges), collect(undirected_edges), n_vertices)
end


@doc raw"""
    label!(G::Graph{T}, edge_labels::Union{Dict{Tuple{Int, Int}, Union{String, Int}}, Nothing}, vertex_labels::Union{Dict{Int, Union{String, Int}}, Nothing}=nothing; name::Symbol=:label) where {T <: Union{Directed, Undirected}}
Given a graph `G`, add labels to the edges and optionally to the vertices with the given `name`.

```jldoctest
julia> G = graph_from_edges(Directed, [[1, 2], [2, 3]])
Directed graph with 3 nodes and the following edges:
(1, 2)(2, 3)

julia> label!(G, Dict((1, 2) => 1), nothing; name=:color)
Directed graph with 3 nodes and the following labeling(s):
label: color
(1, 2) -> 1
(2, 3) -> 0

julia> G = graph_from_labeled_edges(Undirected, Dict((1, 2) => 1, (2, 3) => 2, (1, 3) => 1), nothing; name=:color)
Undirected graph with 3 nodes and the following labeling(s):
label: color
(2, 1) -> 1
(3, 1) -> 1
(3, 2) -> 2

julia> label!(G, nothing, Dict(1 => 1); name=:shading)
Undirected graph with 3 nodes and the following labeling(s):
label: color
(2, 1) -> 1
(3, 1) -> 1
(3, 2) -> 2
label: shading
1 -> 1
2 -> 0
3 -> 0

```
"""
function label!(G::Graph{T},
                    edge_labels::Dict{NTuple{2, Int}, S},
                    vertex_labels::Dict{Int, U};
                    name::Symbol=:label) where {S <: Union{Int, String}, U <: Union{Int, String}, T <: Union{Directed, Undirected}}
  EM = EdgeMap{T, S}(pm_object(G))
  NM = NodeMap{T, U}(pm_object(G))
  set_attribute!(G, name, GraphMap(G, EM, NM))
  for (k, v) in edge_labels
    getproperty(G,name)[k] = v
  end
  for (k, v) in vertex_labels
    @req k <= number_of_vertices(G) "Cannot label a vertex that is not in the graph"
    getproperty(G,name)[k] = v
  end
  return G
end

function label!(G::Graph{T},
                    edge_labels::Dict{NTuple{2, Int}, S},
                    vertex_labels::Nothing;
                    name::Symbol=:label) where {S <: Union{Int, String}, T <: Union{Directed, Undirected}}
  EM = EdgeMap{T, S}(pm_object(G))
  set_attribute!(G, name, GraphMap(G, EM, nothing))
  for (k, v) in edge_labels
    getproperty(G,name)[k] = v
  end
  return G
end

function label!(G::Graph{T},
                    edge_labels::Nothing,
                    vertex_labels::Dict{Int, U};
                    name::Symbol=:label) where {U <: Union{Int, String}, T <: Union{Directed, Undirected}}
  NM = NodeMap{T, U}(pm_object(G))
  set_attribute!(G, name, GraphMap(G, nothing, NM))
  for (k, v) in vertex_labels
    @req k <= number_of_vertices(G) "Cannot label a vertex that is not in the graph"
    getproperty(G,name)[k] = v
  end
  return G
end

@doc raw"""
    graph_from_labeled_edges(edge_labels::Dict{NTuple{Int}, S}, vertex_labels::Union{Nothing, Dict{Int, S}}=nothing; name::Symbol=:label, n_vertices::Int=-1)
    graph_from_labeled_edges(::Type{T}, edge_labels::Dict{NTuple{Int}, S}, vertex_labels::Union{Nothing, Dict{Int, S}}=nothing; name::Symbol=:label, n_vertices::Int=-1) where {T <:Union{Directed, Undirected}, S, U}

Create a graph with a labeling on the edges and optionally vertices.  The graph is constructed from the edges and an optional number of vertices, see [`graph_from_edges`](@ref).
The default name of the labeling on the graph is `label` but any `Symbol` can be passed using the `name` keyword argument.
The labeling can be accessed as a property of the graph, the property is exactly the name passed to the `name` keyword argument.
See [`label!`](@ref) to add additional labels to the graph.

# Examples
```jldoctest
julia> G = graph_from_labeled_edges(Directed, Dict((1, 2) => 1, (2, 3) => 4))
Directed graph with 3 nodes and the following labeling(s):
label: label
(1, 2) -> 1
(2, 3) -> 4

julia> G.label[1, 2]
1

julia> edge = collect(edges(G))[2]
Edge(2, 3)

julia> G.label[edge]
4

julia> K = graph_from_labeled_edges(Dict((1, 2) => "blue", (2, 3) => "green"),
                                    Dict(1 => "red", 2 => "red", 3 => "yellow"); name=:color)
Undirected graph with 3 nodes and the following labeling(s):
label: color
(2, 1) -> blue
(3, 2) -> green
1 -> red
2 -> red
3 -> yellow

julia> K.color[1]
"red"
```
"""
function graph_from_labeled_edges(::Type{T},
                                   edge_labels::Dict{NTuple{2, Int}, <: Union{Int, String}},
                                   vertex_labels::Union{Dict{Int, <: Union{Int, String}}, Nothing}=nothing;
                                   name::Symbol=:label, 
                                   n_vertices::Int=-1) where T <: Union{Directed, Undirected}
  edges = collect(keys(edge_labels))
  G = graph_from_edges(T, edges, n_vertices)
  label!(G, edge_labels, vertex_labels; name=name)
end

function graph_from_labeled_edges(edge_labels::Dict{NTuple{2, Int}, <: Union{Int, String}},
                                   vertex_labels::Union{Dict{Int, <: Union{Int, String}}, Nothing}=nothing;
                                   name::Symbol=:label, n_vertices::Int=-1)
  graph_from_labeled_edges(Undirected, edge_labels, vertex_labels; name=name, n_vertices=n_vertices)
end

@doc raw"""
    adjacency_matrix(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return an unsigned (boolean) adjacency matrix representing a graph `g`. If `g`
is undirected, the adjacency matrix will be symmetric. For `g` being directed,
the adjacency matrix has a `1` at `(u,v)` if there is an edge `u->v`.

# Examples
Adjacency matrix for a directed graph:
```jldoctest
julia> G = Graph{Directed}(3)
Directed graph with 3 nodes and no edges

julia> add_edge!(G,1,3)
true

julia> add_edge!(G,1,2)
true

julia> adjacency_matrix(G)
3×3 IncidenceMatrix
 [2, 3]
 []
 []

julia> matrix(ZZ, adjacency_matrix(G))
[0   1   1]
[0   0   0]
[0   0   0]
```

Adjacency matrix for an undirected graph:
```jldoctest
julia> G = vertex_edge_graph(cube(2))
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(4, 2)(4, 3)

julia> adjacency_matrix(G)
4×4 IncidenceMatrix
 [2, 3]
 [1, 4]
 [1, 4]
 [2, 3]

julia> matrix(ZZ, adjacency_matrix(G))
[0   1   1   0]
[1   0   0   1]
[1   0   0   1]
[0   1   1   0]
```
"""
adjacency_matrix(g::Graph) = Polymake.call_function(:common, Symbol("IncidenceMatrix::new"), nothing, Polymake.common.adjacency_matrix(pm_object(g)))


@doc raw"""
    laplacian_matrix(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the Laplacian matrix of the graph `g`. The Laplacian matrix of a graph
can be written as the difference of `D`, where `D` is a quadratic matrix with
the degrees of `g` on the diagonal, and the adjacency matrix of `g`. For an
undirected graph, the Laplacian matrix is symmetric.

# Examples
```
julia> G = vertex_edge_graph(cube(2))
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(4, 2)(4, 3)

julia> laplacian_matrix(G)
[ 2   -1   -1    0]
[-1    2    0   -1]
[-1    0    2   -1]
[ 0   -1   -1    2]
```
"""
function laplacian_matrix(g::Graph)
  D = diagonal_matrix(degree(g))
  A = matrix(ZZ, adjacency_matrix(g))
  return D-A
end

@doc raw"""
    is_bipartite(g::Graph{Undirected})

Return true if the undirected graph `g` is bipartite.

# Examples
```jldoctest
julia> g = graph_from_edges([[1,2],[2,3],[3,4]]);

julia> is_bipartite(g)
true
```
"""
function is_bipartite(g::Graph{Undirected})
  return Polymake.graph.Graph{Undirected}(ADJACENCY=pm_object(g)).BIPARTITE::Bool
end

@doc raw"""
    maximal_cliques(g::Graph{Undirected})

Returns the maximal cliques of a graph `g` as a `Set{Set{Int}}`.

# Examples
```jldoctest
julia> g = graph_from_edges([[1, 2], [2, 3], [1, 3], [2, 4], [3, 4]])
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 2)(4, 3)

julia> typeof(maximal_cliques(g))
Set{Set{Int64}}

julia> sort.(collect.(maximal_cliques(g)))
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [2, 3, 4]
```
"""
function maximal_cliques(g::Graph{Undirected})
  Set{Set{Int}}(Polymake.to_one_based_indexing.(
    Polymake.call_function(:graph,:max_cliques,g.pm_graph)
  ))
end
