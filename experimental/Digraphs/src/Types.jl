abstract type GAPDigraph <: AbstractGraph{Directed} end

@doc raw"""
    Digraph

A digraph with vertices `1, ..., n`, backed by the GAP package Digraphs.
Every `Digraph` stores the underlying GAP object together with the numbers
of vertices and edges and a flag indicating mutability; see
[`digraph`](@ref) for the constructors, and [`nv`](@ref) and [`ne`](@ref)
for the basic attributes.

```jldoctest
julia> D = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges
```
"""
@attributes mutable struct Digraph <: GAPDigraph
    X::GapObj
    nv::Int64
    ne::Int64
    mut::Bool

    function Digraph(d::GapObj)
        @req DigraphWrap.IsDigraph(d) "the GAP object must be a digraph"
        nv = DigraphWrap.DigraphNrVertices(d)
        ne = DigraphWrap.DigraphNrEdges(d)
        mut = DigraphWrap.IsMutableDigraph(d)
        z = new(d, nv, ne, mut)
        return z
    end
end

GAP.@install GapObj(D::Digraph) = D.X

function Base.show(io::IO, ::MIME"text/plain", D::Digraph)
    @show_name(io, D)
    @show_special(io, D)
    if is_terse(io)
        print(io, "Digraph")
    else
        print(io, "Digraph with ", nv(D), nv(D) == 1 ? " vertex, " : " vertices, ",
              ne(D), ne(D) == 1 ? " edge" : " edges")
    end
end

# ############################################################################
# Conversions between `Digraph` and the Oscar `Graph` types
# (`Graph{T}`, `MixedGraph`) from `src/Combinatorics/Graphs/structs.jl`.
#
# Both directions avoid per-edge crossings of the GAP/Polymake boundary:
# * `Graph -> Digraph` collects the out-neighbours with one Polymake call per
#   vertex and then makes a single GAP `Digraph` call.
# * `Digraph -> Graph` reads the adjacency lists with a single GAP call and
#   constructs the Polymake graph in one call from an `IncidenceMatrix`.
# ############################################################################

# Build a Polymake graph from 1-based adjacency lists in a single call.
function _pm_graph(T::Type{<:Union{Directed, Undirected}}, adj::Vector{Vector{Int}})
    isempty(adj) && return Polymake.Graph{T}(0)
    return Polymake.Graph{T}(Polymake.IncidenceMatrix(adj))
end

@doc raw"""
    digraph(G::Graph{T}; mut::Bool=false) where {T <: Union{Directed, Undirected}}

Return the digraph underlying the Oscar graph `G`. A directed graph is
converted vertex-by-vertex, an undirected graph yields the symmetric digraph
in which every undirected edge `{u, v}` is represented by the two arcs
`u -> v` and `v -> u`. Self-loops are preserved. The result is immutable
unless `mut` is `true`.

# Examples
```jldoctest
julia> G = graph_from_edges(Directed, [[1, 2], [2, 3]])
Directed graph with 3 nodes and the following edges:
(1, 2)(2, 3)

julia> digraph(G)
Digraph with 3 vertices, 2 edges
```
"""
function digraph(G::Graph{T}; mut::Bool=false) where {T <: Union{Directed, Undirected}}
    pm = pm_object(G)
    n = Polymake.nv(pm)
    adj = Vector{Vector{Int}}(undef, n)
    for v in 1:n
        adj[v] = Polymake._outneighbors(pm, v - 1) .+ 1
    end
    filt = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
    return Digraph(DigraphWrap.Digraph(filt, GapObj(adj, recursive=true)))
end

@doc raw"""
    digraph(mg::MixedGraph; mut::Bool=false)

Return the digraph obtained from the mixed graph `mg` by keeping its directed
edges and representing every undirected edge `{u, v}` by the two arcs
`u -> v` and `v -> u`. The result is immutable unless `mut` is `true`.

# Examples
```jldoctest
julia> mg = graph(Mixed, 3);

julia> add_edge!(mg, Directed, 1, 2);

julia> add_edge!(mg, Undirected, 2, 3);

julia> digraph(mg)
Digraph with 3 vertices, 3 edges
```
"""
function digraph(mg::MixedGraph; mut::Bool=false)
    d = Digraph(DigraphWrap.DigraphEdgeUnion(
        GapObj(digraph(_directed_component(mg))),
        GapObj(digraph(_undirected_component(mg)))))
    return mut ? digraph_mutable_copy(d) : digraph_immutable_copy(d)
end

@doc raw"""
    graph(::Type{T}, D::Digraph) where {T <: Union{Directed, Undirected}}

Return the Oscar graph underlying the digraph `D`. For `Directed` the arcs
of `D` are kept; for `Undirected` the digraph must be symmetric and every
pair of opposite arcs becomes one undirected edge. Multiple edges cannot be
represented by `Graph` and raise an error. Self-loops are preserved.

# Examples
```jldoctest
julia> D = digraph([[2, 3], [1, 4], [1], [2]])
Digraph with 4 vertices, 6 edges

julia> graph(Directed, D)
Directed graph with 4 nodes and the following edges:
(1, 2)(1, 3)(2, 1)(2, 4)(3, 1)(4, 2)
```
"""
function graph(::Type{T}, D::Digraph) where {T <: Union{Directed, Undirected}}
    @req !is_multi(D) "the digraph has multiple edges which cannot be represented by an Oscar graph"
    if T === Undirected
        @req is_symmetric(D) "the digraph is not symmetric and cannot be converted to an undirected graph"
    end
    return Graph{T}(_pm_graph(T, out_neighbours(D)))
end

@doc raw"""
    graph(::Type{Mixed}, D::Digraph)

Return the mixed graph obtained from the digraph `D` by moving every pair of
opposite arcs `u -> v`, `v -> u` (and every self-loop) into the undirected
component and keeping all remaining arcs in the directed component. Multiple
edges cannot be represented by `MixedGraph` and raise an error.
"""
function graph(::Type{Mixed}, D::Digraph)
    @req !is_multi(D) "the digraph has multiple edges which cannot be represented by an Oscar graph"
    sym = maximal_symmetric_subdigraph(D)
    adjD = out_neighbours(D)
    adjS = out_neighbours(sym)
    adjR = Vector{Vector{Int}}(undef, length(adjD))
    for i in eachindex(adjD)
        adjR[i] = filter(j -> !(j in adjS[i]), adjD[i])
    end
    gd = Graph{Directed}(_pm_graph(Directed, adjR))
    gu = Graph{Undirected}(_pm_graph(Undirected, out_neighbours(sym)))
    mg = MixedGraph(nv(D))
    mg.directed_component = gd
    mg.undirected_component = gu
    return mg
end

@doc raw"""
    Digraph(G::Graph{T}; mut::Bool=false) where {T <: Union{Directed, Undirected}}

Type-constructor form of [`digraph`](@ref) for Oscar graphs.
"""
Digraph(G::Graph{T}; mut::Bool=false) where {T <: Union{Directed, Undirected}} = digraph(G; mut=mut)

@doc raw"""
    Digraph(mg::MixedGraph; mut::Bool=false)

Type-constructor form of [`digraph`](@ref) for mixed graphs.
"""
Digraph(mg::MixedGraph; mut::Bool=false) = digraph(mg; mut=mut)

@doc raw"""
    Graph{Directed}(D::Digraph)

Type-constructor form of [`graph`](@ref) for directed graphs.
"""
Graph{Directed}(D::Digraph) = graph(Directed, D)

@doc raw"""
    Graph{Undirected}(D::Digraph)

Type-constructor form of [`graph`](@ref) for undirected graphs.
"""
Graph{Undirected}(D::Digraph) = graph(Undirected, D)

@doc raw"""
    MixedGraph(D::Digraph)

Type-constructor form of [`graph`](@ref) for mixed graphs.
"""
MixedGraph(D::Digraph) = graph(Mixed, D)
