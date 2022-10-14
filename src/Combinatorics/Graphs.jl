using Markdown
using JSON
import Oscar: Polyhedron, Polymake, pm_object
import Oscar.Polymake: Directed, Undirected

export
    Directed,
    Edge,
    Graph,
    Undirected,
    add_edge!,
    add_vertex!,
    add_vertices!,
    all_neighbors,
    complete_graph,
    complete_bipartite_graph,
    dst,
    dualgraph,
    edgegraph,
    edges,
    has_edge,
    has_vertex,
    inneighbors,
    ne,
    neighbors,
    nv,
    outneighbors,
    rem_edge!,
    rem_vertex!,
    reverse,
    shortest_path_dijkstra,
    src

################################################################################
################################################################################
##  Constructing and modifying
################################################################################
################################################################################
struct Graph{T <: Union{Directed, Undirected}}
    pm_graph::Polymake.Graph{T}
end
function pm_object(G::Graph{T}) where {T <: Union{Directed, Undirected}}
    return G.pm_graph
end



@doc Markdown.doc"""
    Graph{T}(nverts::Int64) where {T <: Union{Directed, Undirected}}

Construct a graph on `nverts` vertices and no edges. `T` indicates whether the
graph should be `Directed` or `Undirected`.

# Examples
Make a directed graph with 5 vertices and print the number of nodes and edges.
```jldoctest
julia> g = Graph{Directed}(5);

julia> nv(g)
5

julia> ne(g)
0
```
"""
function Graph{T}(nverts::Int64) where {T <: Union{Directed, Undirected}}
    pmg = Polymake.Graph{T}(nverts)
    return Graph{T}(pmg)
end


@doc Markdown.doc"""
    add_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}

Add edge `(s,t)` to the graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> ne(g)
1
```
"""
function add_edge!(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
    Polymake._add_edge(pm_object(g), source-1, target-1)
end


@doc Markdown.doc"""
    rem_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}

Remove edge `(s,t)` from the graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> add_edge!(g, 1, 2);

julia> ne(g)
1

julia> rem_edge!(g, 1, 2);

julia> ne(g)
0
```
"""
function rem_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    return Polymake._rem_edge(pmg, s-1, t-1)
end


@doc Markdown.doc"""
    add_vertex!(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Add a vertex to the graph `g`. The return value is the new vertex.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> nv(g)
2

julia> add_vertex!(g)
3

julia> nv(g)
3
```
"""
function add_vertex!(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    return Polymake._add_vertex(pmg) + 1
end


@doc Markdown.doc"""
    rem_vertex!(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}

Remove the vertex `v` from the graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> nv(g)
2

julia> rem_vertex!(g, 1)

julia> nv(g)
1
```
"""
function rem_vertex!(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    result = Polymake._rem_vertex(pmg, v-1)
    Polymake._squeeze(pmg)
    return result
end


@doc Markdown.doc"""
    add_vertices!(g::Graph{T}, n::Int64) where {T <: Union{Directed, Undirected}}

Add a `n` new vertices to the graph `g`.

# Examples
```jldoctest
julia> g = Graph{Directed}(2);

julia> nv(g)
2

julia> add_vertices!(g, 5);

julia> nv(g)
7
```
"""
function add_vertices!(g::Graph{T}, n::Int64) where {T <: Union{Directed, Undirected}}
    for i = 1:n
        add_vertex!(g)
    end
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


@doc Markdown.doc"""
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


@doc Markdown.doc"""
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


@doc Markdown.doc"""
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


struct EdgeIterator
    pm_itr::Polymake.GraphEdgeIterator{T} where {T <: Union{Directed, Undirected}}
    l::Int64
end
Base.length(eitr::EdgeIterator) = eitr.l
Base.eltype(::Type{EdgeIterator}) = Edge

function Base.iterate(eitr::EdgeIterator, index = 1)
    if index > eitr.l
        return nothing
    else
        e = Polymake.get_element(eitr.pm_itr)
        s = Polymake.first(e)
        t = Polymake.last(e)
        edge = Edge(s+1, t+1)
        Polymake.increment(eitr.pm_itr)
        return (edge, index+1)
    end
end


################################################################################
################################################################################
##  Accessing properties
################################################################################
################################################################################
@doc Markdown.doc"""
    nv(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the number of vertices of a graph.

# Examples
The edge graph of the cube has eight vertices, just like the cube itself.
```jldoctest
julia> c = cube(3);

julia> g = edgegraph(c);

julia> nv(g)
8
```
"""
function nv(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return Polymake.nv(pm_object(g))
end

@doc Markdown.doc"""
    ne(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the number of edges of a graph.

# Examples
The edge graph of the cube has 12 edges just like the cube itself.
```jldoctest
julia> c = cube(3);

julia> g = edgegraph(c);

julia> ne(g)
12
```
"""
function ne(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return Polymake.ne(pm_object(g))
end


@doc Markdown.doc"""
    edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return an iterator over the edges of the graph `g`.

# Examples
A triangle has three edges.
```jldoctest
julia> triangle = simplex(2);

julia> g = edgegraph(triangle);

julia> collect(edges(g))
3-element Vector{Edge}:
 Edge(2, 1)
 Edge(3, 1)
 Edge(3, 2)
```
"""
function edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return EdgeIterator(Polymake.edgeiterator(pm_object(g)), ne(g))
end


@doc Markdown.doc"""
    has_edge(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}

Check for an edge in a graph.

# Examples
Check for the edge $1\to 2$ in the edge graph of a triangle.
```jldoctest
julia> triangle = simplex(2);

julia> g = edgegraph(triangle);

julia> has_edge(g, 1, 2)
true
```
"""
function has_edge(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
    pmg = pm_object(g)
    return Polymake._has_edge(pmg, source-1, target-1)
end
function has_edge(g::Graph{T}, e::Edge) where {T <: Union{Directed, Undirected}}
    return has_edge(g, e.source, e.target)
end


@doc Markdown.doc"""
    has_vertex(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}

Check for a vertex in a graph.

# Examples
The edge graph of a triangle only has 3 vertices.
```jldoctest
julia> triangle = simplex(2);

julia> g = edgegraph(triangle);

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


@doc Markdown.doc"""
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


@doc Markdown.doc"""
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


@doc Markdown.doc"""
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


@doc Markdown.doc"""
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


################################################################################
################################################################################
##  Higher order algorithms
################################################################################
################################################################################
@doc Markdown.doc"""
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
    return _pm_arr_arr_to_group_generators(result)
end


@doc Markdown.doc"""
    automorphism_group(g::Graph{T}) where {T <: Union{Directed, Undirected}}

Return the automorphism group of the graph `g`.

# Examples
```jldoctest
julia> g = complete_graph(4);

julia> automorphism_group(g)
Group([ (3,4), (2,3), (1,2) ])
```
"""
function automorphism_group(g::Graph{T}) where {T <: Union{Directed, Undirected}}
    return _gens_to_group(automorphism_group_generators(g))
end


@doc Markdown.doc"""
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


################################################################################
################################################################################
##  Standard constructions
################################################################################
################################################################################
@doc Markdown.doc"""
    edgegraph(p::Polyhedron)

Return the edge graph of a `Polyhedron`, vertices of the graph correspond to
vertices of the polyhedron, there is an edge between two vertices if the
polyhedron has an edge between the corresponding vertices. The resulting graph
is `Undirected`.

# Examples
Construct the edge graph of the cube. Like the cube it has 8 vertices and 12
edges.
```jldoctest
julia> c = cube(3);

julia> g = edgegraph(c);

julia> nv(g)
8

julia> ne(g)
12
```
"""
function edgegraph(p::Polyhedron)
    pmg = pm_object(p).GRAPH.ADJACENCY
    return Graph{Undirected}(pmg)
end

@doc Markdown.doc"""
    dualgraph(p::Polyhedron)

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

julia> g = dualgraph(c);

julia> nv(g)
6

julia> ne(g)
12
```
"""
function dualgraph(p::Polyhedron)
    pmg = pm_object(p).DUAL_GRAPH.ADJACENCY
    return Graph{Undirected}(pmg)
end



@doc Markdown.doc"""
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


@doc Markdown.doc"""
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



@doc Markdown.doc"""
    visualize(G::Graph{T}) where {T <: Union{Polymake.Directed, Polymake.Undirected}}

Visualize a graph.
"""
function visualize(G::Graph{T}) where {T <: Union{Polymake.Directed, Polymake.Undirected}}
    BigGraph = Polymake.graph.Graph(ADJACENCY=pm_object(G))
    Polymake.visual(BigGraph)
end



# Some standard polytopes from graphs
@doc Markdown.doc"""
    fractional_cut_polytope(G::Graph{Undirected})

Construct the fractional cut polytope of the graph $G$.


# Examples
```jldoctest
julia> G = complete_graph(4);

julia> fractional_cut_polytope(G)
A polyhedron in ambient dimension 6
```
"""
fractional_cut_polytope(G::Graph{Undirected}) = Polyhedron(Polymake.polytope.fractional_cut_polytope(pm_object(G)))


@doc Markdown.doc"""
    fractional_matching_polytope(G::Graph{Undirected})

Construct the fractional matching polytope of the graph $G$.


# Examples
```jldoctest
julia> G = complete_graph(4);

julia> fractional_matching_polytope(G)
A polyhedron in ambient dimension 6
```
"""
fractional_matching_polytope(G::Graph{Undirected}) = Polyhedron(Polymake.polytope.fractional_matching_polytope(pm_object(G)))
