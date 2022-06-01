```@meta
CurrentModule = Oscar.Graphs
```

```@setup oscar
using Oscar
using Oscar.Graphs
```

```@contents
Pages = ["graphs.md"]
```

# Graphs

## Introduction

Graphs are a fundamental object within all of mathematics and computer science.
A *graph* consists of two sets of data:

- a finite set $V := \{1,\ldots,n\}$ of *vertices*; and
- a finite set $E \subseteq V\times V$ of *edges*.

There are two types of graphs, *directed* and *undirected*. For a *directed
graph* the elements of $E$ are considered to be ordered pairs, for an
*undirected graph* the elements of $E$ are unordered pairs or rather sets with
two elements.

The interface is modeled alongside the
[Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/) interface to
allow for easier integration elsewhere.

!!! warning
    The mechanism for removing a vertex is slightly different in out
    implementation to the `Graphs.jl` implementation: In `Graphs.jl` first
    the vertex to be removed is swapped with the last vertex, then the last
    vertex is removed. In our implementation, the vertex is removed and all
    subsequent vertices have their labels changed. Hence edges can be different
    in the two implementations after removing a vertex.

## Construction

```@docs
Graph{T}(nverts::Int64) where {T <: Union{Directed, Undirected}}
dualgraph(p::Polyhedron)
edgegraph(p::Polyhedron)
```

### Modifying graphs
```@docs
add_edge!(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
add_vertices!(g::Graph{T}, n::Int64) where {T <: Union{Directed, Undirected}}
add_vertex!(g::Graph{T}) where {T <: Union{Directed, Undirected}}
rem_edge!(g::Graph{T}, s::Int64, t::Int64) where {T <: Union{Directed, Undirected}}
rem_vertex!(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
```

## Auxiliary functions
```@docs
all_neighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
automorphisms(g::Graph{T}) where {T <: Union{Directed, Undirected}}
complete_graph(n::Int64)
complete_bipartite_graph(n::Int64, m::Int64)
edges(g::Graph{T}) where {T <: Union{Directed, Undirected}}
has_edge(g::Graph{T}, source::Int64, target::Int64) where {T <: Union{Directed, Undirected}}
has_vertex(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
inneighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
ne(g::Graph{T}) where {T <: Union{Directed, Undirected}}
neighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
nv(g::Graph{T}) where {T <: Union{Directed, Undirected}}
outneighbors(g::Graph{T}, v::Int64) where {T <: Union{Directed, Undirected}}
shortest_path_dijkstra
```

### Edges
```@docs
dst(e::Edge)
reverse(e::Edge)
src(e::Edge)
```

## Saving and loading

Objects of type `Graph` can be saved to a file and loaded with the methods
`load` and `save`.  The file is in JSON format and contains the underlying
polymake object. In particular, this file can now be read by both polymake and
Oscar.

