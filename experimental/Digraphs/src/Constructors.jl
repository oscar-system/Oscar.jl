@doc raw"""

    digraph(adj::Vector{Vector{Int}}; mut::Bool=false) -> Digraph
    digraph(labels::Vector{<:AbstractString}, source::Vector{<:AbstractString}, range::Vector{<:AbstractString}; mut::Bool=false) -> Digraph
    digraph(n::Int64, source::Vector{Int}, range::Vector{Int}; mut::Bool=false) -> Digraph
    digraph(list::Vector{<:Any}, func::Function; mut::Bool=false) -> Digraph
    digraph(G::T, list::Vector{<:Any}, act::Function, adj::Function; mut::Bool=false) where T<:GAPGroup -> Digraph
    digraph(obj::GapObj) -> Digraph
    digraph(name::T) where T <: AbstractString -> Digraph

Construct a directed graph (digraph) using one of several convenient interfaces.
By default, an **immutable** digraph is returned; pass `mut=true` to obtain a
**mutable** digraph.

The first form accepts an adjacency list `adj` where the `i`-th entry lists the
out-neighbours of vertex `i`. Repeated entries create multiple edges.

The second form accepts duplicate-free vertex `labels` together with `source`
and `range` vertex-label arrays specifying each edge.

The third form accepts a vertex count `n` together with integer `source` and
`range` arrays for the edges. All entries must be in `[1, n]`.

The fourth form accepts arbitrary objects as vertices and a binary predicate
`func`; an edge from `i` to `j` exists iff `func(list[i], list[j])` returns
`true`.

The fifth form takes a group `G` acting on the vertex list via `act` and an
adjacency predicate `adj` invariant under that action. The action is stored in
the `DigraphGroup` attribute and can accelerate operations such as
`DigraphDiameter`.

The sixth form converts an existing GAP object `obj` (a GAP digraph, a GRAPE
package graph record, or a binary relation on `[1 .. n]`) into a `Digraph`.

The seventh form looks up a named digraph from the built-in database by
`name`. Names are case- and whitespace-insensitive; valid names include
`Diamond`, `Folkman`, and `Brinkmann`. The `graph` suffix may be omitted.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> p = ["a", "b", "c"];

julia> d = digraph(p, p[[1, 1, 2, 2, 3, 3]], p[[2, 3, 1, 3, 1, 2]])
Digraph with 3 vertices, 6 edges

julia> d = digraph(3, [1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2])
Digraph with 3 vertices, 6 edges

julia> rel(x, y) = x != y;

julia> d = digraph([1, 2, 3], rel)
Digraph with 3 vertices, 6 edges

julia> act(x, g) = x^g;

julia> d = digraph(symmetric_group(3), [1, 2, 3], act, rel)
Digraph with 3 vertices, 6 edges

julia> d = digraph("Diamond")
Digraph with 4 vertices, 10 edges

julia> d = GAP.Globals.CompleteGraph(GapObj(symmetric_group(3)))
GAP: rec( adjacencies := [ [ 2, 3 ] ], group := Sym( [ 1 .. 3 ] ), isGraph := true, isSimple := true, order := 3, representatives := [ 1 ], schreierVector := [ -1, 1, 1 ] )

julia> digraph(d)
Digraph with 3 vertices, 6 edges
```
"""
function digraph(adj::Vector{Vector{Int}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(adj, recursive = true))
  return Digraph(d)
end

function digraph(labels::Vector{<:AbstractString}, source::Vector{<:AbstractString}, range::Vector{<:AbstractString}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(labels, recursive = true), GapObj(source, recursive = true), GapObj(range, recursive = true))
  return Digraph(d)
end

function digraph(n::Int64, source::Vector{Int}, range::Vector{Int}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(n), GapObj(source), GapObj(range))
  return Digraph(d)
end

function digraph(list::Vector{<:Any}, func::Function; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(list, recursive = true), GapObj(func))
  return Digraph(d)
end

function digraph(G::T, list::Vector{<:Any}, act::Function, adj::Function; mut::Bool=false) where T<:GAPGroup
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(G), GapObj(list, recursive = true), GapObj(act), GapObj(adj))
  return Digraph(d)
end

function digraph(obj::GapObj)
  d = DigraphWrap.Digraph(obj)
  return Digraph(d)
end

function digraph(name::T) where T <: AbstractString
  d = DigraphWrap.Digraph(GapObj(name))
  return Digraph(d)
end

@doc raw"""
    null_digraph(n::Int64; mut::Bool=false) -> Digraph

Construct a null digraph (a digraph with `n` vertices and no edges).

# Example
```jldoctest
julia> null_digraph(3)
Digraph with 3 vertices, 0 edges
```
"""
function null_digraph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.NullDigraph(filter, n))
end

@doc raw"""
    complete_digraph(n::Int64; mut::Bool=false) -> Digraph

Construct the complete digraph with `n` vertices (every ordered pair of distinct
vertices is an edge, and there are no loops).

# Example
```jldoctest
julia> complete_digraph(3)
Digraph with 3 vertices, 6 edges
```
"""
function complete_digraph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CompleteDigraph(filter, n))
end

@doc raw"""
    complete_bipartite_digraph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct the complete bipartite digraph with parts of sizes `m` and `n`.
Every vertex in one part is connected to every vertex in the other part, and
there are no edges within a part.

# Examples
```jldoctest
julia> complete_bipartite_digraph(2, 3)
Digraph with 5 vertices, 12 edges
```
"""
function complete_bipartite_digraph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CompleteBipartiteDigraph(filter, m, n))
end

@doc raw"""
    cycle_digraph(n::Int64; mut::Bool=false) -> Digraph

Construct the cycle digraph with `n` vertices.

# Example
```jldoctest
julia> cycle_digraph(4)
Digraph with 4 vertices, 4 edges
```
"""
function cycle_digraph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CycleDigraph(filter, n))
end

@doc raw"""
    empty_digraph(n::Int64; mut::Bool=false) -> Digraph

Construct the empty (null) digraph with `n` vertices and no edges.

`NullDigraph` is a synonym for `EmptyDigraph` in GAP; use `null_digraph`
for the equivalent Julia function.

# Example
```jldoctest
julia> empty_digraph(5)
Digraph with 5 vertices, 0 edges
```
"""
function empty_digraph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.EmptyDigraph(filter, n))
end

@doc raw"""
    chain_digraph(n::Int64; mut::Bool=false) -> Digraph

Construct the chain digraph with `n` vertices and `n-1` edges, where there is
a directed edge from `i` to `i+1` for each `i < n`.

# Example
```jldoctest
julia> chain_digraph(4)
Digraph with 4 vertices, 3 edges
```
"""
function chain_digraph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.ChainDigraph(filter, n))
end

@doc raw"""
    johnson_digraph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct the Johnson graph `J(n, k)`, a symmetric digraph whose vertices are
the `k`-subsets of `[1..n]`. Two vertices are adjacent iff their intersection
has size `k-1`.

# Example
```jldoctest
julia> johnson_digraph(4, 2)
Digraph with 6 vertices, 24 edges
```
"""
function johnson_digraph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.JohnsonDigraph(filter, n, k))
end
@doc raw"""
    digraph_from_edges(edges::Vector{Vector{Int64}}; mut::Bool=false) -> Digraph
    digraph_from_edges(edges::Vector{Tuple{Int64,Int64}}; mut::Bool=false) -> Digraph
    digraph_from_edges(n::Int64, edges::Vector{Vector{Int64}}; mut::Bool=false) -> Digraph
    digraph_from_edges(n::Int64, edges::Vector{Tuple{Int64,Int64}}; mut::Bool=false) -> Digraph

Construct a digraph from a list of edges. Each edge is a pair `[source, range]`.
If `n` is given, the digraph will have exactly `n` vertices; otherwise the
vertex count is inferred from the maximum vertex index in `edges`.

# Examples
```jldoctest
julia> digraph_from_edges([[1, 2], [2, 3], [3, 1]])
Digraph with 3 vertices, 3 edges

julia> digraph_from_edges(4, [(1, 2), (2, 3)])
Digraph with 4 vertices, 2 edges
```
"""
function digraph_from_edges(n::Int64, edges::Vector{Tuple{Int64,Int64}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByEdges(filter, GapObj(edges, recursive = true), GapObj(n)))
end

function digraph_from_edges(n::Int64, edges::Vector{Vector{Int64}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByEdges(filter, GapObj(edges, recursive = true), GapObj(n)))
end

function digraph_from_edges(edges::Vector{Tuple{Int64,Int64}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByEdges(filter, GapObj(edges, recursive = true)))
end

function digraph_from_edges(edges::Vector{Vector{Int64}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByEdges(filter, GapObj(edges, recursive = true)))
end

@doc raw"""
    digraph_from_adjacency_matrix(A::Matrix{Int64}; mut::Bool=false) -> Digraph
    digraph_from_adjacency_matrix(A::Matrix{Bool}; mut::Bool=false) -> Digraph

Construct a digraph from its adjacency matrix.

For an integer matrix `A`, entry `A[i, j]` gives the number of edges from
vertex `i` to vertex `j`. For a boolean matrix `A`, entry `A[i, j] == true`
indicates an edge from `i` to `j`.

# Examples
```jldoctest
julia> digraph_from_adjacency_matrix([0 1 0; 1 0 1; 0 1 0])
Digraph with 3 vertices, 4 edges

julia> digraph_from_adjacency_matrix([true false; false true])
Digraph with 2 vertices, 2 edges
```
"""
function digraph_from_adjacency_matrix(A::Matrix{Int64}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByAdjacencyMatrix(filter, GapObj(A)))
end

function digraph_from_adjacency_matrix(A::Matrix{Bool}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByAdjacencyMatrix(filter, GapObj(A)))
end

@doc raw"""
    digraph_from_in_neighbours(inadj::Vector{Vector{Int64}}; mut::Bool=false) -> Digraph
    digraph_from_in_neighbours(inadj::Vector{Tuple{Int64, Int64}}; mut::Bool=false) -> Digraph

Construct a digraph from a list of in-neighbour adjacencies. The `i`-th entry
lists the source vertices of edges whose target is vertex `i`. That is, there
is an edge from `j` to `i` iff `j` is in `inadj[i]`.

# Examples
```jldoctest
julia> digraph_from_in_neighbours([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges
```
"""
function digraph_from_in_neighbours(inadj::Vector{Vector{Int64}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByInNeighbours(filter, GapObj(inadj, recursive = true)))
end

function digraph_from_in_neighbors(inadj::Vector{Tuple{Int64, Int64}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.DigraphByInNeighbours(filter, GapObj(inadj, recursive = true)))
end

@doc raw"""
    random_digraph(n::Int64; mut::Symbol=:mut) -> Digraph
    random_digraph(n::Int64, p::AbstractFloat; mut::Symbol=:mut) -> Digraph

Return a random digraph on `n` vertices, each edge with probability `p`.
Use `mut` for mutability or structure: `:mut`, `:immut`, `:connected`,
`:symmetric`, `:acyclic`, `:eulerian`, or `:hamiltonian`.

# Examples
```jldoctest
julia> random_digraph(100)
Digraph with 100 vertices, 9467 edges

julia> random_digraph(100, 0.5)
Digraph with 100 vertices, 5059 edges

julia> random_digraph(100, 0.25; mut=:connected)
Digraph with 100 vertices, 2552 edges
```
"""
function random_digraph(n::Int64; mut::Symbol=:mut)
  filter = if mut === :mut
    GAP.Globals.IsMutableDigraph
  elseif mut === :immut
    GAP.Globals.IsImmutableDigraph
  elseif mut === :connected
    GAP.Globals.IsConnectedDigraph
  elseif mut === :symmetric
    GAP.Globals.IsSymmetricDigraph
  elseif mut === :acyclic
    GAP.Globals.IsAcyclicDigraph
  elseif mut === :eulerian
    GAP.Globals.IsEulerianDigraph
  elseif mut === :hamiltonian
    GAP.Globals.IsHamiltonianDigraph
  else
    error("unsupported mut $mut")
  end
  return Digraph(DigraphWrap.RandomDigraph(filter, n))
end

function random_digraph(n::Int64, p::AbstractFloat; mut::Symbol=:mut)
  filter = if mut === :mut
    GAP.Globals.IsMutableDigraph
  elseif mut === :immut
    GAP.Globals.IsImmutableDigraph
  elseif mut === :connected
    GAP.Globals.IsConnectedDigraph
  elseif mut === :symmetric
    GAP.Globals.IsSymmetricDigraph
  elseif mut === :acyclic
    GAP.Globals.IsAcyclicDigraph
  elseif mut === :eulerian
    GAP.Globals.IsEulerianDigraph
  elseif mut === :hamiltonian
    GAP.Globals.IsHamiltonianDigraph
  else
    error("unsupported mut $mut")
  end
  r = Rational(p)
  return Digraph(DigraphWrap.RandomDigraph(filter, n, GapObj(p)))
end

@doc raw"""
    random_multi_digraph(n::Int64[, m::Int64]) -> Digraph

Return a random multidigraph with `n` vertices. If `m` is given, the digraph
will have approximately `m` edges; otherwise the number of edges is chosen
randomly.

# Examples
```jldoctest
julia> d = random_multi_digraph(100)
Digraph with 100 vertices, 1124 edges

julia> d = random_multi_digraph(100, 500)
Digraph with 100 vertices, 500 edges
```
"""
function random_multi_digraph(n::Int64)
  return Digraph(DigraphWrap.RandomMultiDigraph(n))
end

function random_multi_digraph(n::Int64, m::Int64)
  return Digraph(DigraphWrap.RandomMultiDigraph(n, m))
end

@doc raw"""
    random_tournament(n::Int64; mut::Bool=false) -> Digraph

Return a random tournament with `n` vertices.

# Examples
```jldoctest
julia> d = random_tournament(10)
Digraph with 10 vertices, 45 edges
```
"""
function random_tournament(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.RandomTournament(filter, n))
end

@doc raw"""
    random_lattice(n::Int64) -> Digraph

Return a random lattice with `m` vertices, where `n <= m <= 2n`.

# Examples
```jldoctest
julia> d = random_lattice(10)
Digraph with 10 vertices, 36 edges
```
"""
function random_lattice(n::Int64)
  return Digraph(DigraphWrap.RandomLattice(n))
end

@doc raw"""
    edge_orbit_digraph(G::T, edges::Vector{Vector{Int64}}; n::Union{Int64,Nothing}=nothing) where T <: GAPGroup -> Digraph
    edge_orbit_digraph(G::T, edges::Vector{Tuple{Int64,Int64}}; n::Union{Int,Nothing}=nothing) where T <: GAPGroup -> Digraph

Construct the edge-orbit digraph of the group `G` acting on the given `edges`.
Each edge is a pair `[source, range]`. If `n` is given, the digraph has exactly
`n` vertices.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> edge_orbit_digraph(G, [[1, 2], [3, 4]])
Digraph with 4 vertices, 12 edges
```
"""
function edge_orbit_digraph(G::T, edges::Vector{Vector{Int64}}; n::Union{Int64,Nothing}=nothing) where T <: GAPGroup
  if n === nothing
    return Digraph(DigraphWrap.EdgeOrbitsDigraph(GapObj(G), GapObj(edges, recursive = true)))
  else
    return Digraph(DigraphWrap.EdgeOrbitsDigraph(GapObj(G), GapObj(edges, recursive = true), n))
  end
end

function edge_orbit_digraph(G::T, edges::Vector{Tuple{Int64,Int64}}; n::Union{Int,Nothing}=nothing) where T <: GAPGroup
  if n === nothing
    return Digraph(DigraphWrap.EdgeOrbitsDigraph(GapObj(G), GapObj(edges, recursive = true)))
  else
    return Digraph(DigraphWrap.EdgeOrbitsDigraph(GapObj(G), GapObj(edges, recursive = true), n))
  end
end

@doc raw"""
    cayley_digraph(G::GAPGroup; gens::Union{Vector{<:GAPGroupElem},Nothing}=nothing; mut::Bool=false) -> Digraph

Construct the Cayley digraph of the group `G` with respect to the generating
set `gens`. If `gens` is not given, the generators of `G` are used.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> cayley_digraph(G)
Digraph with 24 vertices, ...
```
"""
function cayley_digraph(G::T; mut::Bool=false) where T <: GAPGroup
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CayleyDigraph(filter, GapObj(G)))
end

function cayley_digraph(G::T, gens::Vector{<:GAPGroupElem}; mut::Bool=false) where T <: GAPGroup
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CayleyDigraph(filter, GapObj(G), GapObj(gens)))
end

@doc raw"""
    list_named_digraphs(s::AbstractString; level::Int64=2) -> Vector{String}

Search the database of named digraphs for names matching the string `s`.

# Examples
```jldoctest
julia> list_named_digraphs("Dia")
5-element Vector{String}:
 "Diamond"
 ...
```
"""
function list_named_digraphs(s::AbstractString; level::Int64=2)
  return Vector{String}(GAP.Globals.ListNamedDigraphs(GapObj(s), level))
end

@doc raw"""
    andrasfai_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Andrasfai Graph.

# Examples
```jldoctest
julia> d = andrasfai_graph(4)
```
"""
function andrasfai_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.AndrasfaiGraph(filter, n))
end

@doc raw"""
    banana_tree(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Banana Tree.

# Examples
```jldoctest
julia> d = banana_tree(4, 4)
```
"""
function banana_tree(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BananaTree(filter, n, k))
end

@doc raw"""
    binary_tree(n::Int64; mut::Bool=false) -> Digraph

Construct a Binary Tree.

# Examples
```jldoctest
julia> d = binary_tree(4)
```
"""
function binary_tree(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BinaryTree(filter, n))
end

@doc raw"""
    binomial_tree_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Binomial Tree Graph.

# Examples
```jldoctest
julia> d = binomial_tree_graph(4)
```
"""
function binomial_tree_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BinomialTreeGraph(filter, n))
end

@doc raw"""
    bishops_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Bishops Graph.

# Examples
```jldoctest
julia> d = bishops_graph(4, 4)
```
"""
function bishops_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BishopsGraph(filter, m, n))
end

@doc raw"""
    bondy_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Bondy Graph.

# Examples
```jldoctest
julia> d = bondy_graph(4)
```
"""
function bondy_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BondyGraph(filter, n))
end

@doc raw"""
    book_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Book Graph.

# Examples
```jldoctest
julia> d = book_graph(4)
```
"""
function book_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BookGraph(filter, n))
end

@doc raw"""
    burnt_pancake_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Burnt Pancake Graph.

# Examples
```jldoctest
julia> d = burnt_pancake_graph(4)
```
"""
function burnt_pancake_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.BurntPancakeGraph(filter, n))
end

@doc raw"""
    circulant_graph(n::Int64, par::Vector{Int}; mut::Bool=false) -> Digraph

Construct a Circulant Graph.

# Examples
```jldoctest
julia> d = circulant_graph(4, [2, 3])
```
"""
function circulant_graph(n::Int64, par::Vector{Int}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CirculantGraph(filter, n, GapObj(par)))
end

@doc raw"""
    complete_multipartite_digraph(orders::Vector{Int}; mut::Bool=false) -> Digraph

Construct a Complete Multipartite Digraph.

# Examples
```jldoctest
julia> d = complete_multipartite_digraph([2, 3])
```
"""
function complete_multipartite_digraph(orders::Vector{Int}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CompleteMultipartiteDigraph(filter, GapObj(orders)))
end

@doc raw"""
    cycle_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Cycle Graph.

# Examples
```jldoctest
julia> d = cycle_graph(4)
```
"""
function cycle_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.CycleGraph(filter, n))
end

@doc raw"""
    gear_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Gear Graph.

# Examples
```jldoctest
julia> d = gear_graph(4)
```
"""
function gear_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.GearGraph(filter, n))
end

@doc raw"""
    generalised_petersen_graph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Generalised Petersen Graph.

# Examples
```jldoctest
julia> d = generalised_petersen_graph(4, 4)
```
"""
function generalised_petersen_graph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.GeneralisedPetersenGraph(filter, n, k))
end

@doc raw"""
    haar_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Haar Graph.

# Examples
```jldoctest
julia> d = haar_graph(4)
```
"""
function haar_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.HaarGraph(filter, n))
end

@doc raw"""
    halved_cube_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Halved Cube Graph.

# Examples
```jldoctest
julia> d = halved_cube_graph(4)
```
"""
function halved_cube_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.HalvedCubeGraph(filter, n))
end

@doc raw"""
    hanoi_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Hanoi Graph.

# Examples
```jldoctest
julia> d = hanoi_graph(4)
```
"""
function hanoi_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.HanoiGraph(filter, n))
end

@doc raw"""
    helm_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Helm Graph.

# Examples
```jldoctest
julia> d = helm_graph(4)
```
"""
function helm_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.HelmGraph(filter, n))
end

@doc raw"""
    hypercube_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Hypercube Graph.

# Examples
```jldoctest
julia> d = hypercube_graph(4)
```
"""
function hypercube_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.HypercubeGraph(filter, n))
end

@doc raw"""
    keller_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Keller Graph.

# Examples
```jldoctest
julia> d = keller_graph(4)
```
"""
function keller_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.KellerGraph(filter, n))
end

@doc raw"""
    kings_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Kings Graph.

# Examples
```jldoctest
julia> d = kings_graph(4, 4)
```
"""
function kings_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.KingsGraph(filter, m, n))
end

@doc raw"""
    kneser_graph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Kneser Graph.

# Examples
```jldoctest
julia> d = kneser_graph(4, 4)
```
"""
function kneser_graph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.KneserGraph(filter, n, k))
end

@doc raw"""
    knights_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Knights Graph.

# Examples
```jldoctest
julia> d = knights_graph(4, 4)
```
"""
function knights_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.KnightsGraph(filter, m, n))
end

@doc raw"""
    lindgren_sousselier_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Lindgren Sousselier Graph.

# Examples
```jldoctest
julia> d = lindgren_sousselier_graph(4)
```
"""
function lindgren_sousselier_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.LindgrenSousselierGraph(filter, n))
end

@doc raw"""
    lollipop_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Lollipop Graph.

# Examples
```jldoctest
julia> d = lollipop_graph(4, 4)
```
"""
function lollipop_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.LollipopGraph(filter, m, n))
end

@doc raw"""
    mobius_ladder_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Mobius Ladder Graph.

# Examples
```jldoctest
julia> d = mobius_ladder_graph(4)
```
"""
function mobius_ladder_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.MobiusLadderGraph(filter, n))
end

@doc raw"""
    mycielski_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Mycielski Graph.

# Examples
```jldoctest
julia> d = mycielski_graph(4)
```
"""
function mycielski_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.MycielskiGraph(filter, n))
end

@doc raw"""
    odd_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Odd Graph.

# Examples
```jldoctest
julia> d = odd_graph(4)
```
"""
function odd_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.OddGraph(filter, n))
end

@doc raw"""
    pancake_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Pancake Graph.

# Examples
```jldoctest
julia> d = pancake_graph(4)
```
"""
function pancake_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.PancakeGraph(filter, n))
end

@doc raw"""
    path_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Path Graph.

# Examples
```jldoctest
julia> d = path_graph(4)
```
"""
function path_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.PathGraph(filter, n))
end

@doc raw"""
    permutation_star_graph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Permutation Star Graph.

# Examples
```jldoctest
julia> d = permutation_star_graph(4, 4)
```
"""
function permutation_star_graph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.PermutationStarGraph(filter, n, k))
end
@doc raw"""
    prism_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Prism Graph.

# Examples
```jldoctest
julia> d = prism_graph(4)
```
"""
function prism_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.PrismGraph(filter, n))
end

@doc raw"""
    queens_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Queens Graph.

# Examples
```jldoctest
julia> d = queens_graph(4, 4)
```
"""
function queens_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.QueensGraph(filter, m, n))
end

@doc raw"""
    rooks_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Rooks Graph.

# Examples
```jldoctest
julia> d = rooks_graph(4, 4)
```
"""
function rooks_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.RooksGraph(filter, m, n))
end

@doc raw"""
    square_grid_graph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Square Grid Graph.

# Examples
```jldoctest
julia> d = square_grid_graph(4, 4)
```
"""
function square_grid_graph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.SquareGridGraph(filter, n, k))
end

@doc raw"""
    stacked_book_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Stacked Book Graph.

# Examples
```jldoctest
julia> d = stacked_book_graph(4, 4)
```
"""
function stacked_book_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.StackedBookGraph(filter, m, n))
end

@doc raw"""
    stacked_prism_graph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Stacked Prism Graph.

# Examples
```jldoctest
julia> d = stacked_prism_graph(4, 4)
```
"""
function stacked_prism_graph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.StackedPrismGraph(filter, n, k))
end

@doc raw"""
    tadpole_graph(m::Int64, n::Int64; mut::Bool=false) -> Digraph

Construct a Tadpole Graph.

# Examples
```jldoctest
julia> d = tadpole_graph(4, 4)
```
"""
function tadpole_graph(m::Int64, n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.TadpoleGraph(filter, m, n))
end

@doc raw"""
    triangular_grid_graph(n::Int64, k::Int64; mut::Bool=false) -> Digraph

Construct a Triangular Grid Graph.

# Examples
```jldoctest
julia> d = triangular_grid_graph(4, 4)
```
"""
function triangular_grid_graph(n::Int64, k::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.TriangularGridGraph(filter, n, k))
end

@doc raw"""
    walsh_hadamard_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Walsh Hadamard Graph.

# Examples
```jldoctest
julia> d = walsh_hadamard_graph(4)
```
"""
function walsh_hadamard_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.WalshHadamardGraph(filter, n))
end

@doc raw"""
    web_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Web Graph.

# Examples
```jldoctest
julia> d = web_graph(4)
```
"""
function web_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.WebGraph(filter, n))
end

@doc raw"""
    wheel_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Wheel Graph.

# Examples
```jldoctest
julia> d = wheel_graph(4)
```
"""
function wheel_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.WheelGraph(filter, n))
end

@doc raw"""
    windmill_graph(n::Int64, m::Int64; mut::Bool=false) -> Digraph

Construct a Windmill Graph.

# Examples
```jldoctest
julia> d = windmill_graph(4, 4)
```
"""
function windmill_graph(n::Int64, m::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.WindmillGraph(filter, n, m))
end

@doc raw"""
    star_graph(n::Int64; mut::Bool=false) -> Digraph

Construct a Star Graph.

# Examples
```jldoctest
julia> d = star_graph(4)
```
"""
function star_graph(n::Int64; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  return Digraph(DigraphWrap.StarGraph(filter, n))
end

