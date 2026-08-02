```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Constructing Digraphs

Digraphs are created with the constructors in this chapter, which mirror the constructors of Chapter 3 of the GAP Digraphs manual. The basic constructor is [`digraph`](@ref); see also [`Digraph`](@ref) for the type itself.

## Creating digraphs

```@docs
Digraph
is_digraph(d::Digraph)
is_mutable_digraph(d::Digraph)
is_immutable_digraph(d::Digraph)
is_cayley_digraph(d::Digraph)
is_digraph_with_adjacency_function(d::Digraph)
digraph_by_out_neighbours_type()
digraph_family()
digraph(adj::Vector{<:AbstractVector}; mut::Bool=false)
digraph_from_adjacency_matrix(A::Union{Matrix{Int}, Matrix{Bool}}; mut::Bool=false)
digraph_from_edges
edge_orbit_digraph(G::GAPGroup, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}; n::Union{Nothing,Integer}=nothing)
digraph_from_in_neighbours(inadj::Vector{Vector{Int}}; mut::Bool=false)
digraph_from_in_neighbors(inadj::Vector{Vector{Int}}; mut::Bool=false)
cayley_digraph(G::GAPGroup; mut::Bool=false)
list_named_digraphs(s::AbstractString; level::Integer=2)
```

## Changing representations

```@docs
as_binary_relation(d::Digraph)
as_digraph(f::GapObj; mut::Bool=false)
to_grape_graph(d::Digraph)
as_grape_graph(d::Digraph)
as_transformation(d::Digraph)
```

## Conversions with Oscar graphs

Digraphs can be converted to and from the directed, undirected and mixed
`Graph` types of Oscar. The conversions are structural: vertices and edges
are carried over, while labels and other attributes are not. Multiple edges
cannot be represented by Oscar's `Graph` type and raise an error; self-loops
are preserved. An undirected graph is converted to the symmetric digraph
containing both arcs for every edge, and converting a digraph to an
undirected graph requires the digraph to be symmetric. When converting a
digraph to a mixed graph, every pair of opposite arcs becomes one undirected
edge and all remaining arcs stay in the directed component.

```@docs
digraph(G::Graph{T}; mut::Bool=false) where {T <: Union{Directed, Undirected}}
digraph(mg::MixedGraph; mut::Bool=false)
graph(::Type{T}, D::Digraph) where {T <: Union{Directed, Undirected}}
graph(::Type{Mixed}, D::Digraph)
Digraph(G::Graph{T}; mut::Bool=false) where {T <: Union{Directed, Undirected}}
Digraph(mg::MixedGraph; mut::Bool=false)
Graph{Directed}(D::Digraph)
Graph{Undirected}(D::Digraph)
MixedGraph(D::Digraph)
```

## New digraphs from old

```@docs
digraph_immutable_copy(d::Digraph)
digraph_mutable_copy(d::Digraph)
digraph_copy_same_mutability(d::Digraph)
digraph_copy(d::Digraph)
digraph_immutable_copy_if_immutable(d::Digraph)
digraph_immutable_copy_if_mutable(d::Digraph)
digraph_mutable_copy_if_mutable(d::Digraph)
digraph_mutable_copy_if_immutable(d::Digraph)
induced_subdigraph(d::Digraph, verts::Vector{<:Integer})
reduced_digraph(d::Digraph)
reduced_digraph_attr(d::Digraph)
maximal_symmetric_subdigraph(d::Digraph)
maximal_symmetric_subdigraph_attr(d::Digraph)
maximal_symmetric_subdigraph_without_loops(d::Digraph)
maximal_symmetric_subdigraph_without_loops_attr(d::Digraph)
maximal_anti_symmetric_subdigraph(d::Digraph)
maximal_anti_symmetric_subdigraph_attr(d::Digraph)
undirected_spanning_forest(d::Digraph)
undirected_spanning_forest_attr(d::Digraph)
undirected_spanning_tree(d::Digraph)
undirected_spanning_tree_attr(d::Digraph)
digraph_shortest_path_spanning_tree(d::Digraph, root::Integer)
quotient_digraph(d::Digraph, p::Vector{<:AbstractVector{<:Integer}})
digraph_reverse(d::Digraph)
digraph_reverse_attr(d::Digraph)
digraph_dual(d::Digraph)
digraph_dual_attr(d::Digraph)
digraph_symmetric_closure(d::Digraph)
digraph_symmetric_closure_attr(d::Digraph)
digraph_transitive_closure(d::Digraph)
digraph_transitive_closure_attr(d::Digraph)
digraph_reflexive_transitive_closure(d::Digraph)
digraph_reflexive_transitive_closure_attr(d::Digraph)
digraph_transitive_reduction(d::Digraph)
digraph_transitive_reduction_attr(d::Digraph)
digraph_reflexive_transitive_reduction(d::Digraph)
digraph_reflexive_transitive_reduction_attr(d::Digraph)
digraph_add_vertex(d::Digraph)
digraph_add_vertices(d::Digraph, m::Integer)
digraph_add_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
digraph_add_edge_orbit(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
digraph_add_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}})
digraph_remove_vertex(d::Digraph, v::Integer)
digraph_remove_vertices(d::Digraph, verts::Vector{<:Integer})
digraph_remove_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
digraph_remove_edge_orbit(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
digraph_remove_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}})
digraph_remove_loops(d::Digraph)
digraph_remove_loops_attr(d::Digraph)
digraph_remove_all_multiple_edges(d::Digraph)
digraph_remove_all_multiple_edges_attr(d::Digraph)
digraph_contract_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
digraph_reverse_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}})
digraph_reverse_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
digraph_disjoint_union(d1::Digraph, d2::Digraph)
digraph_edge_union(d1::Digraph, d2::Digraph)
digraph_join(d1::Digraph, d2::Digraph)
digraph_cartesian_product(d1::Digraph, d2::Digraph)
digraph_direct_product(d1::Digraph, d2::Digraph)
conormal_product(d1::Digraph, d2::Digraph)
homomorphic_product(d1::Digraph, d2::Digraph)
lexicographic_product(d1::Digraph, d2::Digraph)
digraph_lexicographic_product(d1::Digraph, d2::Digraph)
modular_product(d1::Digraph, d2::Digraph)
digraph_modular_product(d1::Digraph, d2::Digraph)
strong_product(d1::Digraph, d2::Digraph)
digraph_strong_product(d1::Digraph, d2::Digraph)
digraph_cartesian_product_projections(d::Digraph)
digraph_direct_product_projections(d::Digraph)
line_digraph(d::Digraph)
edge_digraph(d::Digraph)
line_undirected_digraph(d::Digraph)
edge_undirected_digraph(d::Digraph)
double_digraph(d::Digraph)
bipartite_double_digraph(d::Digraph)
digraph_add_all_loops(d::Digraph)
digraph_add_all_loops_attr(d::Digraph)
distance_digraph(d::Digraph, i::Integer)
digraph_closure(d::Digraph, k::Integer)
digraph_mycielskian(d::Digraph)
digraph_mycielskian_attr(d::Digraph)
```

## Random digraphs

```@docs
random_digraph(n::Integer; mut::Symbol=:immut)
random_multi_digraph(n::Integer)
random_tournament(n::Integer; mut::Bool=false)
random_lattice(n::Integer; mut::Bool=false)
```

## Standard examples

```@docs
andrasfai_graph(n::Integer; mut::Bool=false)
banana_tree(n::Integer, k::Integer; mut::Bool=false)
binary_tree(m::Integer; mut::Bool=false)
binomial_tree_graph(n::Integer; mut::Bool=false)
bishops_graph(m::Integer, n::Integer; mut::Bool=false)
bishop_graph(m::Integer, n::Integer; mut::Bool=false)
bondy_graph(n::Integer; mut::Bool=false)
book_graph(m::Integer; mut::Bool=false)
burnt_pancake_graph(n::Integer; mut::Bool=false)
pancake_graph(n::Integer; mut::Bool=false)
stacked_book_graph(m::Integer, n::Integer; mut::Bool=false)
chain_digraph(n::Integer; mut::Bool=false)
circulant_graph(n::Integer, par::Vector{<:Integer}; mut::Bool=false)
complete_digraph(n::Integer; mut::Bool=false)
complete_bipartite_digraph(m::Integer, n::Integer; mut::Bool=false)
complete_multipartite_digraph(orders::Vector{<:Integer}; mut::Bool=false)
cycle_digraph(n::Integer; mut::Bool=false)
digraph_cycle(n::Integer; mut::Bool=false)
cycle_symmetric_digraph(n::Integer; mut::Bool=false)
empty_digraph(n::Integer; mut::Bool=false)
null_digraph(n::Integer; mut::Bool=false)
gear_graph(n::Integer; mut::Bool=false)
haar_graph(n::Integer; mut::Bool=false)
halved_cube_graph(n::Integer; mut::Bool=false)
hanoi_graph(n::Integer; mut::Bool=false)
helm_graph(n::Integer; mut::Bool=false)
hypercube_graph(n::Integer; mut::Bool=false)
johnson_digraph(n::Integer, k::Integer; mut::Bool=false)
keller_graph(n::Integer; mut::Bool=false)
kings_graph(m::Integer, n::Integer; mut::Bool=false)
kneser_graph(n::Integer, k::Integer; mut::Bool=false)
knights_graph(m::Integer, n::Integer; mut::Bool=false)
lindgren_sousselier_graph(n::Integer; mut::Bool=false)
lollipop_graph(m::Integer, n::Integer; mut::Bool=false)
mobius_ladder_graph(n::Integer; mut::Bool=false)
mycielski_graph(n::Integer; mut::Bool=false)
odd_graph(n::Integer; mut::Bool=false)
path_graph(n::Integer; mut::Bool=false)
permutation_star_graph(n::Integer, k::Integer; mut::Bool=false)
petersen_digraph(; mut::Bool=false)
generalised_petersen_graph(n::Integer, k::Integer; mut::Bool=false)
prism_graph(n::Integer; mut::Bool=false)
stacked_prism_graph(n::Integer, k::Integer; mut::Bool=false)
queens_graph(m::Integer, n::Integer; mut::Bool=false)
queen_graph(m::Integer, n::Integer; mut::Bool=false)
rooks_graph(m::Integer, n::Integer; mut::Bool=false)
rook_graph(m::Integer, n::Integer; mut::Bool=false)
square_grid_graph(n::Integer, k::Integer; mut::Bool=false)
grid_graph(n::Integer, k::Integer; mut::Bool=false)
triangular_grid_graph(n::Integer, k::Integer; mut::Bool=false)
star_graph(k::Integer; mut::Bool=false)
tadpole_graph(m::Integer, n::Integer; mut::Bool=false)
walsh_hadamard_graph(n::Integer; mut::Bool=false)
web_graph(n::Integer; mut::Bool=false)
wheel_graph(n::Integer; mut::Bool=false)
windmill_graph(n::Integer, m::Integer; mut::Bool=false)
```
