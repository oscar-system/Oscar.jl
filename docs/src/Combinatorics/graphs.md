```@meta
CurrentModule = Oscar
```

# Graphs

## Constructors

```@docs
Graph(::Type{T}, n::Int) where {T<:Union{Directed, Undirected}}
Graph(::Type{T}, A::MatElem) where {T<:Union{Directed, Undirected}}
random_graph(n::Int, p::AbstractFloat, T::Type)
empty_graph(::Type{T}, n::Int) where {T<:Union{Directed, Undirected}}
complete_graph(n::Int64)
complete_bipartite_graph(n::Int64, m::Int64)
petersen_graph()
clebsch_graph()
cayley_graph
cayley_graph_index_map
cayley_graph_vertex
```            

### Others

```@docs
cartesian_product(g1::Graph, g2::Graph)
complement(g::Graph)
line_graph(g::Graph)
induced_subgraph(g::Graph, v::Vector{Int})
union(g1::Graph{T}, g2::Graph{T}) where {T<:Union{Directed, Undirected}}
```

## Properties

```@docs
n_vertices(g::Graph)
n_edges(g::Graph)
vertices(g::Graph)
edges(g::Graph)
is_directed(g::Graph{Directed})
is_undirected(g::Graph{Undirected})
neighbors(g::Graph{T}, v::Int) where {T<:Union{Directed, Undirected}}
in_degree(g::Graph{T}, v::Int) where {T<:Union{Directed, Undirected}}
out_degree(g::Graph{T}, v::Int) where {T<:Union{Directed, Undirected}}
degree(g::Graph{T}, v::Int) where {T<:Union{Directed, Undirected}}
is_adjacent(g::Graph, u::Int, v::Int)
common_neighbors(g::Graph, u::Int, v::Int)
```

## Modifying the graph

```@docs
add_vertex!(g::Graph)
add_vertices!(g::Graph, n::Int)
add_edge!(g::Graph, u::Int, v::Int)
add_edges!(g::Graph, edges::Vector)
remove_vertex!(g::Graph, v::Int)
remove_edge!(g::Graph, u::Int, v::Int)
```

## Algorithms

```@docs
connected_components(g::Graph{Undirected})
is_connected(g::Graph{Undirected})
is_bipartite(g::Graph)
bipartition(g::Graph)
dijkstra_shortest_path(g::Graph, source::Int)
is_acyclic(g::Graph{Directed})
topological_sort(g::Graph{Directed})
is_strongly_connected(g::Graph{Directed})
strongly_connected_components(g::Graph{Directed})
girth(g::Graph)
diameter(g::Graph)
radius(g::Graph)
center(g::Graph)
is_tree(g::Graph)
spanning_tree(g::Graph{Undirected})
a_star(g::Graph, source::Int, target::Int)
eulerian_trail(g::Graph{Undirected})
eulerian_circuit(g::Graph{Undirected})
maximal_cliques(g::Graph{Undirected})
labelings(G::Graph)
has_disjoint_automorphisms(G::Graph)
disjoint_automorphisms(G::Graph)
```

### Edges

```@docs
src(e::Edge)
dst(e::Edge)
reverse(e::Edge)
```

### Auto-generated from Oscar

```@docs
is_regular(g::Graph)
is_complete(g::Graph)
is_biconnected(g::Graph)
laplacian_matrix(g::Graph)
is_vertex_transitive(g::Graph)
is_edge_transitive(g::Graph)
is_arc_transitive(g::Graph)
is_distance_transitive(g::Graph)
is_distance_regular(g::Graph)
```
