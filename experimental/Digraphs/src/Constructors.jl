# ============================================================================
# Constructors.jl - full wrapper for the GAP Digraphs package, Chapter 3
# Reference: https://docs.gap-system.org/pkg/digraphs/doc/chap3_mj.html
#
# This file follows the exact order of the GAP manual chapter "Creating
# digraphs" (3.1 -> 3.2 -> 3.3 -> 3.4 -> 3.5) and wraps every overload of
# every operation, including the `*Attr` attribute variants and the official
# GAP aliases.
#
# Conventions:
# * The keyword argument `mut::Bool=false` is only provided for GAP
#   operations whose signature contains the optional filter `[filt, ]`; a
#   `false` value selects `IsImmutableDigraph` and `true` selects
#   `IsMutableDigraph`.
# * All other operations (mostly the in-place operations of section 3.3)
#   follow the GAP semantics: the mutability of the result is determined by
#   the input digraph (a mutable digraph is changed in-place, an immutable
#   digraph is left unchanged and a new digraph is returned).
# ============================================================================

# Helper: return the appropriate GAP filter based on mutability.
function _filt(mut::Bool)
    return mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
end

# Helpers for converting Julia edge data to GAP lists.
_gap_edge(e::Tuple{Int,Int}) = GapObj([e[1], e[2]])
_gap_edge(e::Vector{Int}) = GapObj(e)
_gap_edges(edges::Vector{Tuple{Int,Int}}) = GapObj([collect(Int, e) for e in edges], recursive=true)
_gap_edges(edges::Vector{Vector{Int}}) = GapObj(edges, recursive=true)
_gap_verts(verts::Vector{<:Integer}) = GapObj(collect(Int, verts))

# ############################################################################
# Section 3.1  Creating digraphs
# ############################################################################

# ----------------------------------------------------------------------------
# 3.1-1 .. 3.1-5  categories
# ----------------------------------------------------------------------------
"""
    is_digraph(d::Digraph) -> Bool

Return `true` if `d` belongs to the GAP category `IsDigraph`. Every digraph
in the GAP Digraphs package does.
"""
function is_digraph(d::Digraph)
    return DigraphWrap.IsDigraph(GapObj(d))
end

"""
    is_mutable_digraph(d::Digraph) -> Bool

Return `true` if `d` is a mutable digraph (i.e. it lies in the GAP category
`IsMutableDigraph`). A mutable digraph may be changed in-place by the
operations of section 3.3.
"""
function is_mutable_digraph(d::Digraph)
    return DigraphWrap.IsMutableDigraph(GapObj(d))
end

"""
    is_immutable_digraph(d::Digraph) -> Bool

Return `true` if `d` is an immutable digraph (i.e. it lies in the GAP category
`IsImmutableDigraph`). Immutable digraphs are attribute-storing; mutating
operations return a new digraph instead of changing `d`.
"""
function is_immutable_digraph(d::Digraph)
    return DigraphWrap.IsImmutableDigraph(GapObj(d))
end

"""
    is_cayley_digraph(d::Digraph) -> Bool

Return `true` if `d` is a Cayley digraph of a group (constructed with
[`cayley_digraph`](@ref)).
"""
function is_cayley_digraph(d::Digraph)
    return DigraphWrap.IsCayleyDigraph(GapObj(d))
end

"""
    is_digraph_with_adjacency_function(d::Digraph) -> Bool

Return `true` if `d` was created using an adjacency function (see the fourth
and fifth variants of [`digraph`](@ref)).
"""
function is_digraph_with_adjacency_function(d::Digraph)
    return DigraphWrap.IsDigraphWithAdjacencyFunction(GapObj(d))
end

# ----------------------------------------------------------------------------
# 3.1-6  DigraphByOutNeighboursType / DigraphFamily
# ----------------------------------------------------------------------------
"""
    digraph_by_out_neighbours_type() -> GapObj

Return the GAP global `DigraphByOutNeighboursType`, the type of all digraphs.
"""
function digraph_by_out_neighbours_type()
    return GAP.Globals.DigraphByOutNeighboursType
end

"""
    digraph_family() -> GapObj

Return the GAP family `DigraphFamily`, the family of all digraphs.
"""
function digraph_family()
    return GAP.Globals.DigraphFamily
end

# ----------------------------------------------------------------------------
# 3.1-7  Digraph  (core constructor, multiple dispatch)
# ----------------------------------------------------------------------------
@doc raw"""
    digraph(adj::Vector{<:AbstractVector}; mut::Bool=false) -> Digraph
    digraph(labels::Vector, source::Vector, range::Vector; mut::Bool=false) -> Digraph
    digraph(n::Integer, source::Vector{<:Integer}, range::Vector{<:Integer}; mut::Bool=false) -> Digraph
    digraph(list::Vector, func::Function; mut::Bool=false) -> Digraph
    digraph(G::GAPGroup, list::Vector, act::Function, adj::Function; mut::Bool=false) -> Digraph
    digraph(obj::GapObj) -> Digraph
    digraph(name::AbstractString) -> Digraph

Construct a digraph using one of several convenient interfaces. By default an
immutable digraph is returned; pass `mut=true` to obtain a mutable digraph
(where the GAP operation accepts a filter).

The first form accepts an adjacency list `adj`: the `i`-th entry lists the
out-neighbours of vertex `i` (positive integers), and an entry occurring `k`
times gives `k` parallel edges. The entries are validated by GAP.

The second form accepts a duplicate-free list `labels` of vertex labels
together with equal-length `source` and `range` lists; the resulting digraph
has vertices `[1 .. length(labels)]` labelled by the elements of `labels`.

The third form accepts a vertex count `n` together with integer `source` and
`range` lists; every entry must lie in `[1 .. n]`.

The fourth form takes a list of objects and a binary predicate `func`; there
is an edge from `i` to `j` iff `func(list[i], list[j])` returns `true`.

The fifth form takes a group `G` acting on `list` via `act` and an adjacency
predicate `adj` which is invariant under the action; the action is stored in
the attribute `DigraphGroup` of the result and speeds up operations such as
the diameter.

The sixth form converts an existing GAP object `obj` (a GAP digraph, a GRAPE
package graph, or a binary relation on `[1 .. n]`) into a `Digraph`; its
mutability is determined by `obj`.

The seventh form looks up a named digraph from the built-in GAP database by
`name` (names are case- and whitespace-insensitive). Named digraphs are
always immutable.

# Examples
```jldoctest
julia> d = digraph([[2, 5, 8, 10], [2, 3, 4, 2, 5, 6, 8, 9, 10], [1],
       [3, 5, 7, 8, 10], [2, 5, 7], [3, 6, 7, 9, 10], [1, 4],
       [1, 5, 9], [1, 2, 7, 8], [3, 5]])
Digraph with 10 vertices, 38 edges

julia> digraph(["a", "b", "c"], ["a"], ["b"])
Digraph with 3 vertices, 1 edges

julia> digraph(5, [1, 2, 2, 4, 1, 1], [2, 3, 5, 5, 1, 1])
Digraph with 5 vertices, 6 edges

julia> digraph(collect(1:10), (x, y) -> true)
Digraph with 10 vertices, 100 edges

julia> digraph("Diamond")
Digraph with 4 vertices, 10 edges
```
"""
function digraph(adj::Vector{<:AbstractVector}; mut::Bool=false)
    d = DigraphWrap.Digraph(_filt(mut), GapObj(adj, recursive=true))
    return Digraph(d)
end

function digraph(labels::Vector, source::Vector, range::Vector; mut::Bool=false)
    @req length(source) == length(range) "source and range must have equal length"
    d = DigraphWrap.Digraph(_filt(mut), GapObj(labels, recursive=true),
                            GapObj(source, recursive=true),
                            GapObj(range, recursive=true))
    return Digraph(d)
end

function digraph(n::Integer, source::Vector{<:Integer}, range::Vector{<:Integer}; mut::Bool=false)
    @req length(source) == length(range) "source and range must have equal length"
    d = DigraphWrap.Digraph(_filt(mut), Int64(n), _gap_verts(source), _gap_verts(range))
    return Digraph(d)
end

function digraph(list::Vector, func::Function; mut::Bool=false)
    d = DigraphWrap.Digraph(_filt(mut), GapObj(list, recursive=true), GapObj(func))
    return Digraph(d)
end

function digraph(G::GAPGroup, list::Vector, act::Function, adj::Function; mut::Bool=false)
    d = DigraphWrap.Digraph(_filt(mut), GapObj(G), GapObj(list, recursive=true),
                            GapObj(act), GapObj(adj))
    return Digraph(d)
end

function digraph(obj::GapObj)
    return Digraph(DigraphWrap.Digraph(obj))
end

function digraph(name::AbstractString)
    return Digraph(DigraphWrap.Digraph(GapObj(String(name))))
end

# ----------------------------------------------------------------------------
# 3.1-8  DigraphByAdjacencyMatrix
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_from_adjacency_matrix(A::Union{Matrix{Int}, Matrix{Bool}}; mut::Bool=false) -> Digraph

Construct a digraph from its adjacency matrix `A`. For an integer matrix,
entry `A[i, j]` is the number of edges from vertex `i` to vertex `j`; for a
boolean matrix, `A[i, j] == true` indicates an edge from `i` to `j`.

# Examples
```jldoctest
julia> digraph_from_adjacency_matrix([0 1 0; 1 0 1; 0 1 0])
Digraph with 3 vertices, 4 edges

julia> digraph_from_adjacency_matrix([true false; false true])
Digraph with 2 vertices, 2 edges
```
"""
function digraph_from_adjacency_matrix(A::Union{Matrix{Int}, Matrix{Bool}}; mut::Bool=false)
    return Digraph(DigraphWrap.DigraphByAdjacencyMatrix(_filt(mut), GapObj(A)))
end

# ----------------------------------------------------------------------------
# 3.1-9  DigraphByEdges
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_from_edges(edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}; mut::Bool=false) -> Digraph
    digraph_from_edges(edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, n::Integer; mut::Bool=false) -> Digraph

Construct a digraph from the list `edges` of `[source, range]` pairs. If the
vertex count `n` is given, the digraph has exactly `n` vertices; otherwise
the smallest possible number of vertices is used.

# Examples
```jldoctest
julia> digraph_from_edges([[1, 2], [2, 3], [3, 1]])
Digraph with 3 vertices, 3 edges

julia> digraph_from_edges([(1, 2), (2, 3), (3, 1)], 5)
Digraph with 5 vertices, 3 edges
```
"""
function digraph_from_edges(edges::Vector{Vector{Int}}; mut::Bool=false)
    return Digraph(DigraphWrap.DigraphByEdges(_filt(mut), GapObj(edges, recursive=true)))
end

function digraph_from_edges(edges::Vector{Vector{Int}}, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.DigraphByEdges(_filt(mut), GapObj(edges, recursive=true), Int64(n)))
end

function digraph_from_edges(edges::Vector{Tuple{Int,Int}}; mut::Bool=false)
    return digraph_from_edges([collect(Int, e) for e in edges]; mut=mut)
end

function digraph_from_edges(edges::Vector{Tuple{Int,Int}}, n::Integer; mut::Bool=false)
    return digraph_from_edges([collect(Int, e) for e in edges], n; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.1-10  EdgeOrbitsDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    edge_orbit_digraph(G::GAPGroup, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}; n::Union{Nothing,Integer}=nothing) -> Digraph

Return the immutable digraph with `n` vertices whose edges are the union of
the orbits of the edges in `edges` under the permutation group `G`. If `n` is
not given, the largest moved point of `G` is used. The result is always
immutable.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> edge_orbit_digraph(G, [[1, 2], [3, 4]])
Digraph with 4 vertices, 12 edges
```
"""
function edge_orbit_digraph(G::GAPGroup, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}; n::Union{Nothing,Integer}=nothing)
    eg = _gap_edges(edges)
    if n === nothing
        return Digraph(DigraphWrap.EdgeOrbitsDigraph(GapObj(G), eg))
    else
        return Digraph(DigraphWrap.EdgeOrbitsDigraph(GapObj(G), eg, Int(n)))
    end
end

# ----------------------------------------------------------------------------
# 3.1-11  DigraphByInNeighbours
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_from_in_neighbours(inadj::Vector{Vector{Int}}; mut::Bool=false) -> Digraph

Construct a digraph from a list of in-neighbour adjacencies: the `i`-th entry
`inadj[i]` lists the sources of all edges whose range is vertex `i`. An entry
occurring `k` times in `inadj[i]` gives `k` parallel edges from that source
to `i`.

# Examples
```jldoctest
julia> digraph_from_in_neighbours([[2], [1, 3], [2]])
Digraph with 3 vertices, 4 edges
```
"""
function digraph_from_in_neighbours(inadj::Vector{Vector{Int}}; mut::Bool=false)
    return Digraph(DigraphWrap.DigraphByInNeighbours(_filt(mut), GapObj(inadj, recursive=true)))
end

"""
    digraph_from_in_neighbors(inadj::Vector{Vector{Int}}; mut::Bool=false) -> Digraph

American-spelling alias for [`digraph_from_in_neighbours`](@ref).
"""
function digraph_from_in_neighbors(inadj::Vector{Vector{Int}}; mut::Bool=false)
    return digraph_from_in_neighbours(inadj; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.1-12  CayleyDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    cayley_digraph(G::GAPGroup; mut::Bool=false) -> Digraph
    cayley_digraph(G::GAPGroup, gens::Vector{<:GAPGroupElem}; mut::Bool=false) -> Digraph

Construct the Cayley digraph of the group `G` with respect to the generating
set `gens`; if `gens` is not given, the generators of `G` are used. The
vertices correspond to the elements of `G` in the order of `Set(G)`, and the
group and generators can be recovered with `group_of_cayley_digraph` and
`generators_of_cayley_digraph`.

# Examples
```jldoctest
julia> cayley_digraph(symmetric_group(4))
Digraph with 24 vertices, 48 edges
```
"""
function cayley_digraph(G::GAPGroup; mut::Bool=false)
    return Digraph(DigraphWrap.CayleyDigraph(_filt(mut), GapObj(G)))
end

function cayley_digraph(G::GAPGroup, gens::Vector{<:GAPGroupElem}; mut::Bool=false)
    return Digraph(DigraphWrap.CayleyDigraph(_filt(mut), GapObj(G), GapObj(gens, recursive=true)))
end

# ----------------------------------------------------------------------------
# 3.1-13  ListNamedDigraphs
# ----------------------------------------------------------------------------
@doc raw"""
    list_named_digraphs(s::AbstractString; level::Integer=2) -> Vector{String}

Search the GAP database of named digraphs for names matching the string `s`.
The `level` controls the search: `level=1` matches only names starting with
`s`, `level=2` (the default) matches names containing `s` as a substring, and
`level=3` performs a substring search ignoring non-alphanumeric characters.
The search is case- and whitespace-insensitive.

# Examples
```jldoctest
julia> length(list_named_digraphs("diamond"))
1
```
"""
function list_named_digraphs(s::AbstractString; level::Integer=2)
    return Vector{String}(DigraphWrap.ListNamedDigraphs(GapObj(String(s)), Int(level)))
end

# ############################################################################
# Section 3.2  Changing representations
# ############################################################################

# ----------------------------------------------------------------------------
# 3.2-1  AsBinaryRelation
# ----------------------------------------------------------------------------
"""
    as_binary_relation(d::Digraph) -> GapObj

If `d` is a digraph with no multiple edges and a positive number `n` of
vertices, return the GAP binary relation on `[1 .. n]` whose pairs are the
edges of `d`.
"""
function as_binary_relation(d::Digraph)
    return DigraphWrap.AsBinaryRelation(GapObj(d))
end

# ----------------------------------------------------------------------------
# 3.2-2  AsDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    as_digraph(f::GapObj; mut::Bool=false) -> Digraph
    as_digraph(f::GapObj, n::Integer; mut::Bool=false) -> Digraph

Convert a GAP transformation, permutation, partial perm, or binary relation
`f` into a digraph on `[1 .. n]` vertices with an edge `v -> v^f` for every
vertex `v`. If `n` is omitted, the degree of the transformation (or the
largest moved point of the permutation, etc.) is used.

An `ArgumentError` is thrown if `v^f` exceeds `n` for some `v`, in which case
GAP returns `fail`.

# Examples
```jldoctest
julia> f = GAP.Globals.Transformation([4, 3, 3, 1, 7, 9, 10, 4, 2, 3]);

julia> as_digraph(f)
Digraph with 10 vertices, 10 edges
```
"""
function as_digraph(f::GapObj; mut::Bool=false)
    result = DigraphWrap.AsDigraph(_filt(mut), f)
    @req result != GAP.Globals.fail "as_digraph: GAP returns fail because the image of some point exceeds the number of vertices; see AsDigraph (Digraphs manual 3.2-2)."
    return Digraph(result)
end

function as_digraph(f::GapObj, n::Integer; mut::Bool=false)
    result = DigraphWrap.AsDigraph(_filt(mut), f, Int64(n))
    @req result != GAP.Globals.fail "as_digraph: GAP returns fail because the image of some point exceeds the number of vertices; see AsDigraph (Digraphs manual 3.2-2)."
    return Digraph(result)
end

# ----------------------------------------------------------------------------
# 3.2-3  Graph  (to GRAPE)
# ----------------------------------------------------------------------------
"""
    to_grape_graph(d::Digraph) -> GapObj

Convert the digraph `d` to a GRAPE package graph isomorphic to `d`. Multiple
edges of a multidigraph are reduced to a single edge, since GRAPE does not
support multiple edges.
"""
function to_grape_graph(d::Digraph)
    return DigraphWrap.Graph(GapObj(d))
end

# ----------------------------------------------------------------------------
# 3.2-4  AsGraph  (to GRAPE, cached)
# ----------------------------------------------------------------------------
"""
    as_grape_graph(d::Digraph) -> GapObj

Same as [`to_grape_graph`](@ref), but for an immutable digraph the result is
stored as an attribute of `d`, so subsequent calls return the same GAP
object.
"""
function as_grape_graph(d::Digraph)
    return DigraphWrap.AsGraph(GapObj(d))
end

# ----------------------------------------------------------------------------
# 3.2-5  AsTransformation
# ----------------------------------------------------------------------------
@doc raw"""
    as_transformation(d::Digraph) -> GapObj

If `d` is a functional digraph (every vertex has a unique out-neighbour),
return the GAP transformation whose image of vertex `i` is that unique
out-neighbour. An `ArgumentError` is thrown otherwise, since GAP returns
`fail` for a non-functional digraph.

# Examples
```jldoctest
julia> as_transformation(digraph([[1], [3], [2]]))
GAP: Transformation( [ 1, 3, 2 ] )
```
"""
function as_transformation(d::Digraph)
    result = DigraphWrap.AsTransformation(GapObj(d))
    @req result != GAP.Globals.fail "as_transformation: GAP returns fail because d is not a functional digraph."
    return result
end

# ############################################################################
# Section 3.3  New digraphs from old
# ############################################################################

# ----------------------------------------------------------------------------
# 3.3-1  DigraphImmutableCopy / DigraphMutableCopy / DigraphCopySameMutability
# ----------------------------------------------------------------------------
"""
    digraph_immutable_copy(d::Digraph) -> Digraph

Return a new immutable copy of `d`, retaining none of its attributes or
properties.
"""
function digraph_immutable_copy(d::Digraph)
    return Digraph(DigraphWrap.DigraphImmutableCopy(GapObj(d)))
end

"""
    digraph_mutable_copy(d::Digraph) -> Digraph

Return a new mutable copy of `d`, retaining none of its attributes or
properties.
"""
function digraph_mutable_copy(d::Digraph)
    return Digraph(DigraphWrap.DigraphMutableCopy(GapObj(d)))
end

"""
    digraph_copy_same_mutability(d::Digraph) -> Digraph

Return a new copy of `d` with the same mutability as `d`, retaining none of
its attributes or properties.
"""
function digraph_copy_same_mutability(d::Digraph)
    return Digraph(DigraphWrap.DigraphCopySameMutability(GapObj(d)))
end

"""
    digraph_copy(d::Digraph) -> Digraph

Synonym for [`digraph_copy_same_mutability`](@ref), mirroring the GAP alias
`DigraphCopy`.
"""
function digraph_copy(d::Digraph)
    return digraph_copy_same_mutability(d)
end

# ----------------------------------------------------------------------------
# 3.3-2  Digraph*CopyIf* (conditional copies)
# ----------------------------------------------------------------------------
"""
    digraph_immutable_copy_if_immutable(d::Digraph) -> Digraph

Return `d` itself if `d` is mutable, and a new immutable copy of `d`
otherwise.
"""
function digraph_immutable_copy_if_immutable(d::Digraph)
    return Digraph(DigraphWrap.DigraphImmutableCopyIfImmutable(GapObj(d)))
end

"""
    digraph_immutable_copy_if_mutable(d::Digraph) -> Digraph

Return a new immutable copy of `d` if `d` is mutable, and `d` itself
otherwise.
"""
function digraph_immutable_copy_if_mutable(d::Digraph)
    return Digraph(DigraphWrap.DigraphImmutableCopyIfMutable(GapObj(d)))
end

"""
    digraph_mutable_copy_if_mutable(d::Digraph) -> Digraph

Return a new mutable copy of `d` if `d` is mutable, and `d` itself otherwise.
"""
function digraph_mutable_copy_if_mutable(d::Digraph)
    return Digraph(DigraphWrap.DigraphMutableCopyIfMutable(GapObj(d)))
end

"""
    digraph_mutable_copy_if_immutable(d::Digraph) -> Digraph

Return a new mutable copy of `d` if `d` is immutable, and `d` itself
otherwise.
"""
function digraph_mutable_copy_if_immutable(d::Digraph)
    return Digraph(DigraphWrap.DigraphMutableCopyIfImmutable(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-3  InducedSubdigraph
# ----------------------------------------------------------------------------
@doc raw"""
    induced_subdigraph(d::Digraph, verts::Vector{<:Integer}) -> Digraph

Return the subdigraph of `d` induced by the vertices `verts`: it retains
precisely the vertices in `verts` and the edges whose source and range are
both in `verts`. If `d` is mutable it is changed in-place; otherwise a new
immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[1, 1, 2, 3, 4, 4], [1, 3, 4], [3, 1], [1, 1]])
Digraph with 4 vertices, 13 edges

julia> induced_subdigraph(d, [1, 3, 4])
Digraph with 3 vertices, 9 edges
```
"""
function induced_subdigraph(d::Digraph, verts::Vector{<:Integer})
    return Digraph(DigraphWrap.InducedSubdigraph(GapObj(d), _gap_verts(verts)))
end

# ----------------------------------------------------------------------------
# 3.3-4  ReducedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    reduced_digraph(d::Digraph) -> Digraph

Return the subdigraph of `d` induced by its non-isolated vertices (the
vertices which are the source or range of some edge). The ordering and the
labels of the remaining vertices are preserved. If `d` is mutable it is
changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[1, 2], [], [], [1, 4], []])
Digraph with 5 vertices, 4 edges

julia> reduced_digraph(d)
Digraph with 3 vertices, 4 edges
```
"""
function reduced_digraph(d::Digraph)
    return Digraph(DigraphWrap.ReducedDigraph(GapObj(d)))
end

"""
    reduced_digraph_attr(d::Digraph) -> Digraph

Same as [`reduced_digraph`](@ref), but for an immutable digraph the result is
stored as an attribute of `d`.
"""
function reduced_digraph_attr(d::Digraph)
    return Digraph(DigraphWrap.ReducedDigraphAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-5  MaximalSymmetricSubdigraph
# ----------------------------------------------------------------------------
@doc raw"""
    maximal_symmetric_subdigraph(d::Digraph) -> Digraph

Return the symmetric digraph on the vertices of `d` obtained by ignoring the
multiplicity of edges and keeping an edge `[u, v]` only if `[v, u]` is also
an edge of `d`. If `d` is mutable it is changed in-place; otherwise a new
immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1, 3], [4], [3, 1]])
Digraph with 4 vertices, 7 edges

julia> maximal_symmetric_subdigraph(d)
Digraph with 4 vertices, 4 edges
```
"""
function maximal_symmetric_subdigraph(d::Digraph)
    return Digraph(DigraphWrap.MaximalSymmetricSubdigraph(GapObj(d)))
end

"""
    maximal_symmetric_subdigraph_attr(d::Digraph) -> Digraph

Same as [`maximal_symmetric_subdigraph`](@ref), but for an immutable digraph
the result is stored as an attribute of `d`.
"""
function maximal_symmetric_subdigraph_attr(d::Digraph)
    return Digraph(DigraphWrap.MaximalSymmetricSubdigraphAttr(GapObj(d)))
end

@doc raw"""
    maximal_symmetric_subdigraph_without_loops(d::Digraph) -> Digraph

Same as [`maximal_symmetric_subdigraph`](@ref), except that loops are removed
from the result.

# Examples
```jldoctest
julia> d = digraph([[1, 1, 2], [1], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> maximal_symmetric_subdigraph_without_loops(d)
Digraph with 3 vertices, 2 edges
```
"""
function maximal_symmetric_subdigraph_without_loops(d::Digraph)
    return Digraph(DigraphWrap.MaximalSymmetricSubdigraphWithoutLoops(GapObj(d)))
end

"""
    maximal_symmetric_subdigraph_without_loops_attr(d::Digraph) -> Digraph

Same as [`maximal_symmetric_subdigraph_without_loops`](@ref), but for an
immutable digraph the result is stored as an attribute of `d`.
"""
function maximal_symmetric_subdigraph_without_loops_attr(d::Digraph)
    return Digraph(DigraphWrap.MaximalSymmetricSubdigraphWithoutLoopsAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-6  MaximalAntiSymmetricSubdigraph
# ----------------------------------------------------------------------------
@doc raw"""
    maximal_anti_symmetric_subdigraph(d::Digraph) -> Digraph

Return the anti-symmetric subdigraph of `d` obtained by discarding duplicate
edges and, for every pair of opposite edges `[i, j]` and `[j, i]` with
distinct vertices, discarding the edge `[max(i, j), min(i, j)]`. If `d` is
mutable it is changed in-place; otherwise a new immutable digraph is
returned.

# Examples
```jldoctest
julia> d = digraph([[2, 2], [1, 3], [4], [3, 1]])
Digraph with 4 vertices, 7 edges

julia> maximal_anti_symmetric_subdigraph(d)
Digraph with 4 vertices, 4 edges
```
"""
function maximal_anti_symmetric_subdigraph(d::Digraph)
    return Digraph(DigraphWrap.MaximalAntiSymmetricSubdigraph(GapObj(d)))
end

"""
    maximal_anti_symmetric_subdigraph_attr(d::Digraph) -> Digraph

Same as [`maximal_anti_symmetric_subdigraph`](@ref), but for an immutable
digraph the result is stored as an attribute of `d`.
"""
function maximal_anti_symmetric_subdigraph_attr(d::Digraph)
    return Digraph(DigraphWrap.MaximalAntiSymmetricSubdigraphAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-7  UndirectedSpanningForest / UndirectedSpanningTree
# ----------------------------------------------------------------------------
@doc raw"""
    undirected_spanning_forest(d::Digraph) -> Digraph

Return an undirected spanning forest of the underlying undirected graph of
`d`. GAP returns `fail` for a digraph without vertices, in which case an
`ArgumentError` is thrown. If `d` is mutable it is changed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 1, 3], [1], [4], [3, 4, 3]])
Digraph with 4 vertices, 9 edges

julia> undirected_spanning_forest(d)
Digraph with 4 vertices, 4 edges
```
"""
function undirected_spanning_forest(d::Digraph)
    result = DigraphWrap.UndirectedSpanningForest(GapObj(d))
    @req result != GAP.Globals.fail "undirected_spanning_forest: GAP returns fail because the digraph has no vertices."
    return Digraph(result)
end

"""
    undirected_spanning_forest_attr(d::Digraph) -> Digraph

Same as [`undirected_spanning_forest`](@ref), but for an immutable digraph
the result is stored as an attribute of `d`.
"""
function undirected_spanning_forest_attr(d::Digraph)
    result = DigraphWrap.UndirectedSpanningForestAttr(GapObj(d))
    @req result != GAP.Globals.fail "undirected_spanning_forest_attr: GAP returns fail because the digraph has no vertices."
    return Digraph(result)
end

@doc raw"""
    undirected_spanning_tree(d::Digraph) -> Digraph

Return an undirected spanning tree of the underlying undirected graph of `d`.
GAP returns `fail` when the digraph has no vertices or when its maximal
symmetric subdigraph is not connected, in which case an `ArgumentError` is
thrown. If `d` is mutable it is changed in-place; otherwise a new immutable
digraph is returned.

# Examples
```jldoctest
julia> d = complete_digraph(4)
Digraph with 4 vertices, 12 edges

julia> undirected_spanning_tree(d)
Digraph with 4 vertices, 6 edges
```
"""
function undirected_spanning_tree(d::Digraph)
    result = DigraphWrap.UndirectedSpanningTree(GapObj(d))
    @req result != GAP.Globals.fail "undirected_spanning_tree: GAP returns fail because the digraph has no vertices or its maximal symmetric subdigraph is not connected."
    return Digraph(result)
end

"""
    undirected_spanning_tree_attr(d::Digraph) -> Digraph

Same as [`undirected_spanning_tree`](@ref), but for an immutable digraph the
result is stored as an attribute of `d`.
"""
function undirected_spanning_tree_attr(d::Digraph)
    result = DigraphWrap.UndirectedSpanningTreeAttr(GapObj(d))
    @req result != GAP.Globals.fail "undirected_spanning_tree_attr: GAP returns fail because the digraph has no vertices or its maximal symmetric subdigraph is not connected."
    return Digraph(result)
end

# ----------------------------------------------------------------------------
# 3.3-8  DigraphShortestPathSpanningTree
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_shortest_path_spanning_tree(d::Digraph, root::Integer) -> Digraph

Return the shortest-path spanning tree of `d` rooted at the vertex `root`.
GAP returns `fail` when some vertex of `d` is not reachable from `root`, in
which case an `ArgumentError` is thrown. If `d` is mutable it is changed
in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = complete_digraph(4)
Digraph with 4 vertices, 12 edges

julia> digraph_shortest_path_spanning_tree(d, 1)
Digraph with 4 vertices, 3 edges
```
"""
function digraph_shortest_path_spanning_tree(d::Digraph, root::Integer)
    result = DigraphWrap.DigraphShortestPathSpanningTree(GapObj(d), Int(root))
    @req result != GAP.Globals.fail "digraph_shortest_path_spanning_tree: GAP returns fail because not all vertices are reachable from the root."
    return Digraph(result)
end

# ----------------------------------------------------------------------------
# 3.3-9  QuotientDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    quotient_digraph(d::Digraph, p::Vector{<:AbstractVector{<:Integer}}) -> Digraph

Return the quotient digraph of `d` with respect to the partition `p` of its
vertices: the vertices of each part are amalgamated into a single vertex. The
quotient has no multiple edges. If `d` is mutable it is changed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2, 1], [4], [1], [1, 3, 4]])
Digraph with 4 vertices, 7 edges

julia> quotient_digraph(d, [[1], [2, 4], [3]])
Digraph with 3 vertices, 6 edges
```
"""
function quotient_digraph(d::Digraph, p::Vector{<:AbstractVector{<:Integer}})
    gap_p = GapObj([collect(Int, part) for part in p], recursive=true)
    return Digraph(DigraphWrap.QuotientDigraph(GapObj(d), gap_p))
end

# ----------------------------------------------------------------------------
# 3.3-10  DigraphReverse
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_reverse(d::Digraph) -> Digraph

Return the reverse of `d`, i.e. the digraph obtained by reversing the
orientation of every edge. If `d` is mutable it is changed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[3], [1, 3, 5], [1], [1, 2, 4], [2, 3, 5]])
Digraph with 5 vertices, 11 edges

julia> digraph_reverse(d)
Digraph with 5 vertices, 11 edges
```
"""
function digraph_reverse(d::Digraph)
    return Digraph(DigraphWrap.DigraphReverse(GapObj(d)))
end

"""
    digraph_reverse_attr(d::Digraph) -> Digraph

Same as [`digraph_reverse`](@ref), but for an immutable digraph the result is
stored as an attribute of `d`.
"""
function digraph_reverse_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphReverseAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-11  DigraphDual
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_dual(d::Digraph) -> Digraph

Return the dual (complement) of `d`: the digraph on the same vertices with an
edge from `i` to `j` iff there is no edge from `i` to `j` in `d`. If `d` is
mutable it is changed in-place; otherwise a new immutable digraph is
returned.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [], [4, 6], [5], [], [7, 8, 9], [], [], []])
Digraph with 9 vertices, 8 edges

julia> digraph_dual(d)
Digraph with 9 vertices, 73 edges
```
"""
function digraph_dual(d::Digraph)
    return Digraph(DigraphWrap.DigraphDual(GapObj(d)))
end

"""
    digraph_dual_attr(d::Digraph) -> Digraph

Same as [`digraph_dual`](@ref), but for an immutable digraph the result is
stored as an attribute of `d`.
"""
function digraph_dual_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphDualAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-12  DigraphSymmetricClosure
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_symmetric_closure(d::Digraph) -> Digraph

Return the smallest symmetric digraph with the same vertices as `d` which
contains all the edges of `d`. If `d` is mutable it is changed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3], [2, 4], [1], [3, 4]])
Digraph with 4 vertices, 8 edges

julia> digraph_symmetric_closure(d)
Digraph with 4 vertices, 11 edges
```
"""
function digraph_symmetric_closure(d::Digraph)
    return Digraph(DigraphWrap.DigraphSymmetricClosure(GapObj(d)))
end

"""
    digraph_symmetric_closure_attr(d::Digraph) -> Digraph

Same as [`digraph_symmetric_closure`](@ref), but for an immutable digraph the
result is stored as an attribute of `d`.
"""
function digraph_symmetric_closure_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphSymmetricClosureAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-13  DigraphTransitiveClosure
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_transitive_closure(d::Digraph) -> Digraph

Return the transitive closure of the digraph `d` (which must have no multiple
edges): the least transitive digraph containing `d`. If `d` is mutable it is
changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_transitive_closure(d)
Digraph with 3 vertices, 3 edges
```
"""
function digraph_transitive_closure(d::Digraph)
    return Digraph(DigraphWrap.DigraphTransitiveClosure(GapObj(d)))
end

"""
    digraph_transitive_closure_attr(d::Digraph) -> Digraph

Same as [`digraph_transitive_closure`](@ref), but for an immutable digraph
the result is stored as an attribute of `d`.
"""
function digraph_transitive_closure_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphTransitiveClosureAttr(GapObj(d)))
end

@doc raw"""
    digraph_reflexive_transitive_closure(d::Digraph) -> Digraph

Return the reflexive transitive closure of the digraph `d` (which must have
no multiple edges): the least reflexive and transitive digraph containing
`d`. If `d` is mutable it is changed in-place; otherwise a new immutable
digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_reflexive_transitive_closure(d)
Digraph with 3 vertices, 6 edges
```
"""
function digraph_reflexive_transitive_closure(d::Digraph)
    return Digraph(DigraphWrap.DigraphReflexiveTransitiveClosure(GapObj(d)))
end

"""
    digraph_reflexive_transitive_closure_attr(d::Digraph) -> Digraph

Same as [`digraph_reflexive_transitive_closure`](@ref), but for an immutable
digraph the result is stored as an attribute of `d`.
"""
function digraph_reflexive_transitive_closure_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphReflexiveTransitiveClosureAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-14  DigraphTransitiveReduction
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_transitive_reduction(d::Digraph) -> Digraph

Return the transitive reduction of the topologically sortable digraph `d`
(which must have no multiple edges): the least subdigraph of `d` whose
transitive closure equals the transitive closure of `d`. If `d` is mutable it
is changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> digraph_transitive_reduction(d)
Digraph with 3 vertices, 2 edges
```
"""
function digraph_transitive_reduction(d::Digraph)
    return Digraph(DigraphWrap.DigraphTransitiveReduction(GapObj(d)))
end

"""
    digraph_transitive_reduction_attr(d::Digraph) -> Digraph

Same as [`digraph_transitive_reduction`](@ref), but for an immutable digraph
the result is stored as an attribute of `d`.
"""
function digraph_transitive_reduction_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphTransitiveReductionAttr(GapObj(d)))
end

@doc raw"""
    digraph_reflexive_transitive_reduction(d::Digraph) -> Digraph

Return the reflexive transitive reduction of the topologically sortable
digraph `d` (which must have no multiple edges): the least subdigraph of `d`
whose reflexive transitive closure equals the reflexive transitive closure of
`d`. If `d` is mutable it is changed in-place; otherwise a new immutable
digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2, 3], [3], []])
Digraph with 3 vertices, 3 edges

julia> digraph_reflexive_transitive_reduction(d)
Digraph with 3 vertices, 2 edges
```
"""
function digraph_reflexive_transitive_reduction(d::Digraph)
    return Digraph(DigraphWrap.DigraphReflexiveTransitiveReduction(GapObj(d)))
end

"""
    digraph_reflexive_transitive_reduction_attr(d::Digraph) -> Digraph

Same as [`digraph_reflexive_transitive_reduction`](@ref), but for an
immutable digraph the result is stored as an attribute of `d`.
"""
function digraph_reflexive_transitive_reduction_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphReflexiveTransitiveReductionAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-15  DigraphAddVertex
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_add_vertex(d::Digraph) -> Digraph
    digraph_add_vertex(d::Digraph, label) -> Digraph

Return the digraph obtained from `d` by adding a single new vertex (and no
edges), optionally labelled by `label`. If `d` is mutable the vertex is added
in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_add_vertex(d)
Digraph with 4 vertices, 6 edges
```
"""
function digraph_add_vertex(d::Digraph)
    return Digraph(DigraphWrap.DigraphAddVertex(GapObj(d)))
end

function digraph_add_vertex(d::Digraph, label)
    return Digraph(DigraphWrap.DigraphAddVertex(GapObj(d), GapObj(label)))
end

# ----------------------------------------------------------------------------
# 3.3-16  DigraphAddVertices
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_add_vertices(d::Digraph, m::Integer) -> Digraph
    digraph_add_vertices(d::Digraph, labels::Vector) -> Digraph

Return the digraph obtained from `d` by adding `m` new vertices, or by adding
one new vertex for each entry of `labels` (labelled accordingly). If `d` is
mutable the vertices are added in-place; otherwise a new immutable digraph is
returned.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_add_vertices(d, 3)
Digraph with 6 vertices, 6 edges
```
"""
function digraph_add_vertices(d::Digraph, m::Integer)
    return Digraph(DigraphWrap.DigraphAddVertices(GapObj(d), Int(m)))
end

function digraph_add_vertices(d::Digraph, labels::Vector)
    return Digraph(DigraphWrap.DigraphAddVertices(GapObj(d), GapObj(labels, recursive=true)))
end

# ----------------------------------------------------------------------------
# 3.3-17  DigraphAddEdge
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_add_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}}) -> Digraph
    digraph_add_edge(d::Digraph, src::Integer, ran::Integer) -> Digraph

Return the digraph obtained from `d` by adding the edge `[src, ran]`. If `d`
is mutable the edge is added in-place; otherwise a new immutable digraph is
returned.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_add_edge(d, [3, 1])
Digraph with 3 vertices, 3 edges
```
"""
function digraph_add_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
    return Digraph(DigraphWrap.DigraphAddEdge(GapObj(d), _gap_edge(edge)))
end

function digraph_add_edge(d::Digraph, src::Integer, ran::Integer)
    return Digraph(DigraphWrap.DigraphAddEdge(GapObj(d), Int(src), Int(ran)))
end

# ----------------------------------------------------------------------------
# 3.3-18  DigraphAddEdgeOrbit
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_add_edge_orbit(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}}) -> Digraph

Return a new digraph with the same vertices and edges as the immutable
digraph `d` together with the orbit of `edge` under the action of the
`DigraphGroup` of `d`. If `edge` is already an edge of `d`, then `d` itself
is returned unchanged.
"""
function digraph_add_edge_orbit(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
    return Digraph(DigraphWrap.DigraphAddEdgeOrbit(GapObj(d), _gap_edge(edge)))
end

# ----------------------------------------------------------------------------
# 3.3-19  DigraphAddEdges
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_add_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}) -> Digraph

Return the digraph obtained from `d` by adding all edges in `edges`; an edge
occurring `k` times is added `k` times. If `d` is mutable the edges are added
in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = empty_digraph(3)
Digraph with 3 vertices, 0 edges

julia> digraph_add_edges(d, [[1, 2], [2, 3]])
Digraph with 3 vertices, 2 edges
```
"""
function digraph_add_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}})
    return Digraph(DigraphWrap.DigraphAddEdges(GapObj(d), _gap_edges(edges)))
end

# ----------------------------------------------------------------------------
# 3.3-20  DigraphRemoveVertex
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_vertex(d::Digraph, v::Integer) -> Digraph

Return the digraph obtained from `d` by removing the vertex `v` together with
all edges incident to `v`. If `d` is mutable the vertex is removed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph(["a", "b", "c"], ["a", "a", "b", "c", "c"], ["b", "c", "a", "a", "c"])
Digraph with 3 vertices, 5 edges

julia> digraph_remove_vertex(d, 2)
Digraph with 2 vertices, 3 edges
```
"""
function digraph_remove_vertex(d::Digraph, v::Integer)
    return Digraph(DigraphWrap.DigraphRemoveVertex(GapObj(d), Int(v)))
end

# ----------------------------------------------------------------------------
# 3.3-21  DigraphRemoveVertices
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_vertices(d::Digraph, verts::Vector{<:Integer}) -> Digraph

Return the digraph obtained from `d` by removing every vertex in `verts`
together with all incident edges. If `d` is mutable the vertices are removed
in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[3], [1, 3, 5], [1], [1, 2, 4], [2, 3, 5]])
Digraph with 5 vertices, 11 edges

julia> digraph_remove_vertices(d, [2, 4])
Digraph with 3 vertices, 4 edges
```
"""
function digraph_remove_vertices(d::Digraph, verts::Vector{<:Integer})
    return Digraph(DigraphWrap.DigraphRemoveVertices(GapObj(d), _gap_verts(verts)))
end

# ----------------------------------------------------------------------------
# 3.3-22  DigraphRemoveEdge
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}}) -> Digraph
    digraph_remove_edge(d::Digraph, src::Integer, ran::Integer) -> Digraph

Return the digraph obtained from `d` by removing the edge `[src, ran]`; `d`
must have no multiple edges. If `d` is mutable the edge is removed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = cycle_digraph(4)
Digraph with 4 vertices, 4 edges

julia> digraph_remove_edge(d, [4, 1])
Digraph with 4 vertices, 3 edges
```
"""
function digraph_remove_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
    return Digraph(DigraphWrap.DigraphRemoveEdge(GapObj(d), _gap_edge(edge)))
end

function digraph_remove_edge(d::Digraph, src::Integer, ran::Integer)
    return Digraph(DigraphWrap.DigraphRemoveEdge(GapObj(d), Int(src), Int(ran)))
end

# ----------------------------------------------------------------------------
# 3.3-23  DigraphRemoveEdgeOrbit
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_edge_orbit(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}}) -> Digraph

Return a new digraph with the same vertices as the immutable digraph `d` and
with the orbit of `edge` under the action of the `DigraphGroup` of `d`
removed. If `edge` is not an edge of `d`, then `d` itself is returned
unchanged.
"""
function digraph_remove_edge_orbit(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
    return Digraph(DigraphWrap.DigraphRemoveEdgeOrbit(GapObj(d), _gap_edge(edge)))
end

# ----------------------------------------------------------------------------
# 3.3-24  DigraphRemoveEdges
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}) -> Digraph

Return the digraph obtained from `d` by removing all edges in `edges`; `d`
must have no multiple edges (unless `edges` is empty). Invalid entries of
`edges` are ignored. If `d` is mutable the edges are removed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = complete_digraph(3)
Digraph with 3 vertices, 6 edges

julia> digraph_remove_edges(d, [[1, 2], [2, 1]])
Digraph with 3 vertices, 4 edges
```
"""
function digraph_remove_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}})
    return Digraph(DigraphWrap.DigraphRemoveEdges(GapObj(d), _gap_edges(edges)))
end

# ----------------------------------------------------------------------------
# 3.3-25  DigraphRemoveLoops
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_loops(d::Digraph) -> Digraph

Return the digraph obtained from `d` by removing every loop (an edge with
equal source and range). If `d` is mutable the loops are removed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 4], [1, 4], [3, 4], [1, 4, 5], [1, 5]])
Digraph with 5 vertices, 12 edges

julia> digraph_remove_loops(d)
Digraph with 5 vertices, 8 edges
```
"""
function digraph_remove_loops(d::Digraph)
    return Digraph(DigraphWrap.DigraphRemoveLoops(GapObj(d)))
end

"""
    digraph_remove_loops_attr(d::Digraph) -> Digraph

Same as [`digraph_remove_loops`](@ref), but for an immutable digraph the
result is stored as an attribute of `d`.
"""
function digraph_remove_loops_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphRemoveLoopsAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-26  DigraphRemoveAllMultipleEdges
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_remove_all_multiple_edges(d::Digraph) -> Digraph

Return the largest subdigraph of `d` with no multiple edges. If `d` is
mutable the multiple edges are removed in-place; otherwise a new immutable
digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[1, 2, 3, 2], [1, 1, 3], [2, 2, 2]])
Digraph with 3 vertices, 10 edges

julia> digraph_remove_all_multiple_edges(d)
Digraph with 3 vertices, 6 edges
```
"""
function digraph_remove_all_multiple_edges(d::Digraph)
    return Digraph(DigraphWrap.DigraphRemoveAllMultipleEdges(GapObj(d)))
end

"""
    digraph_remove_all_multiple_edges_attr(d::Digraph) -> Digraph

Same as [`digraph_remove_all_multiple_edges`](@ref), but for an immutable
digraph the result is stored as an attribute of `d`.
"""
function digraph_remove_all_multiple_edges_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphRemoveAllMultipleEdgesAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-27  DigraphContractEdge
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_contract_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}}) -> Digraph
    digraph_contract_edge(d::Digraph, src::Integer, ran::Integer) -> Digraph

Return the digraph obtained from the digraph `d` (which must have no multiple
edges) by contracting the edge `[src, ran]` with `src != ran`: the two
vertices are merged, incident edges keep their direction, and vertex labels
of the merged vertices are combined into a list. If `d` is mutable it is
changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph_from_edges([[1, 2], [2, 1]])
Digraph with 2 vertices, 2 edges

julia> digraph_contract_edge(d, 1, 2)
Digraph with 1 vertices, 0 edges
```
"""
function digraph_contract_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
    return Digraph(DigraphWrap.DigraphContractEdge(GapObj(d), _gap_edge(edge)))
end

function digraph_contract_edge(d::Digraph, src::Integer, ran::Integer)
    return Digraph(DigraphWrap.DigraphContractEdge(GapObj(d), Int(src), Int(ran)))
end

# ----------------------------------------------------------------------------
# 3.3-28  DigraphReverseEdges / DigraphReverseEdge
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_reverse_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}) -> Digraph

Return the digraph obtained from the digraph `d` (which must have no multiple
edges) by reversing the orientation of every edge in `edges`. The result may
have multiple edges. If `d` is mutable the edges are reversed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_reverse_edges(d, [[1, 2]])
Digraph with 3 vertices, 2 edges
```
"""
function digraph_reverse_edges(d::Digraph, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}})
    return Digraph(DigraphWrap.DigraphReverseEdges(GapObj(d), _gap_edges(edges)))
end

@doc raw"""
    digraph_reverse_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}}) -> Digraph
    digraph_reverse_edge(d::Digraph, src::Integer, ran::Integer) -> Digraph

Reverse a single edge of the digraph `d` (which must have no multiple edges);
see [`digraph_reverse_edges`](@ref).

# Examples
```jldoctest
julia> d = digraph([[2], [3], []])
Digraph with 3 vertices, 2 edges

julia> digraph_reverse_edge(d, 2, 3)
Digraph with 3 vertices, 2 edges
```
"""
function digraph_reverse_edge(d::Digraph, edge::Union{Vector{Int}, Tuple{Int,Int}})
    return Digraph(DigraphWrap.DigraphReverseEdge(GapObj(d), _gap_edge(edge)))
end

function digraph_reverse_edge(d::Digraph, src::Integer, ran::Integer)
    return Digraph(DigraphWrap.DigraphReverseEdge(GapObj(d), Int(src), Int(ran)))
end

# ----------------------------------------------------------------------------
# 3.3-29  DigraphDisjointUnion
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_disjoint_union(d1::Digraph, d2::Digraph) -> Digraph
    digraph_disjoint_union(d1::Digraph, d2::Digraph, rest::Digraph...) -> Digraph
    digraph_disjoint_union(ds::Vector{Digraph}) -> Digraph

Return the disjoint union of the given digraphs: the vertex set and the edge
list are the disjoint unions of those of the arguments, with the vertices of
the `i`-th digraph shifted by the total number of vertices of the preceding
ones. Previously set vertex labels are lost. If the first digraph is mutable
it is changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> digraph_disjoint_union(cycle_digraph(3), cycle_digraph(4))
Digraph with 7 vertices, 7 edges

julia> digraph_disjoint_union([cycle_digraph(3), cycle_digraph(4)])
Digraph with 7 vertices, 7 edges
```
"""
function digraph_disjoint_union(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.DigraphDisjointUnion(GapObj(d1), GapObj(d2)))
end

function digraph_disjoint_union(ds::Vector{Digraph})
    @req !isempty(ds) "digraph_disjoint_union requires a non-empty list of digraphs"
    gap_list = GapObj([GapObj(d) for d in ds], recursive=true)
    return Digraph(DigraphWrap.DigraphDisjointUnion(gap_list))
end

function digraph_disjoint_union(d1::Digraph, d2::Digraph, rest::Digraph...)
    return digraph_disjoint_union(Digraph[d1, d2, rest...])
end

# ----------------------------------------------------------------------------
# 3.3-30  DigraphEdgeUnion
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_edge_union(d1::Digraph, d2::Digraph) -> Digraph
    digraph_edge_union(d1::Digraph, d2::Digraph, rest::Digraph...) -> Digraph
    digraph_edge_union(ds::Vector{Digraph}) -> Digraph

Return the edge union of the given digraphs: the vertex set is the union of
the vertex sets (its size is the maximum over the arguments) and the edge
list is the concatenation of the edge lists. Previously set vertex labels are
lost. If the first digraph is mutable it is changed in-place; otherwise a new
immutable digraph is returned.

# Examples
```jldoctest
julia> d1 = digraph([[2], [1]])
Digraph with 2 vertices, 2 edges

julia> d2 = digraph([[2, 3], [2], [1]])
Digraph with 3 vertices, 4 edges

julia> digraph_edge_union(d1, d2)
Digraph with 3 vertices, 6 edges
```
"""
function digraph_edge_union(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.DigraphEdgeUnion(GapObj(d1), GapObj(d2)))
end

function digraph_edge_union(ds::Vector{Digraph})
    @req !isempty(ds) "digraph_edge_union requires a non-empty list of digraphs"
    gap_list = GapObj([GapObj(d) for d in ds], recursive=true)
    return Digraph(DigraphWrap.DigraphEdgeUnion(gap_list))
end

function digraph_edge_union(d1::Digraph, d2::Digraph, rest::Digraph...)
    return digraph_edge_union(Digraph[d1, d2, rest...])
end

# ----------------------------------------------------------------------------
# 3.3-31  DigraphJoin
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_join(d1::Digraph, d2::Digraph) -> Digraph
    digraph_join(d1::Digraph, d2::Digraph, rest::Digraph...) -> Digraph
    digraph_join(ds::Vector{Digraph}) -> Digraph

Return the join of the given digraphs: the disjoint union together with every
edge between vertices belonging to different arguments. Previously set vertex
labels are lost. If the first digraph is mutable it is changed in-place;
otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> digraph_join(complete_digraph(3), cycle_digraph(3))
Digraph with 6 vertices, 27 edges
```
"""
function digraph_join(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.DigraphJoin(GapObj(d1), GapObj(d2)))
end

function digraph_join(ds::Vector{Digraph})
    @req !isempty(ds) "digraph_join requires a non-empty list of digraphs"
    gap_list = GapObj([GapObj(d) for d in ds], recursive=true)
    return Digraph(DigraphWrap.DigraphJoin(gap_list))
end

function digraph_join(d1::Digraph, d2::Digraph, rest::Digraph...)
    return digraph_join(Digraph[d1, d2, rest...])
end

# ----------------------------------------------------------------------------
# 3.3-32  DigraphCartesianProduct
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_cartesian_product(d1::Digraph, d2::Digraph) -> Digraph
    digraph_cartesian_product(d1::Digraph, d2::Digraph, rest::Digraph...) -> Digraph
    digraph_cartesian_product(ds::Vector{Digraph}) -> Digraph

Return the Cartesian product of the given digraphs. The vertices of the
result are encoded as integers, and the original pairs are recovered as
vertex labels: `digraph_vertex_label(d, i)` is the pair of labels of the two
factors. The projections onto the coordinates are available via
[`digraph_cartesian_product_projections`](@ref).

# Examples
```jldoctest
julia> digraph_cartesian_product(chain_digraph(4), cycle_digraph(3))
Digraph with 12 vertices, 21 edges
```
"""
function digraph_cartesian_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.DigraphCartesianProduct(GapObj(d1), GapObj(d2)))
end

function digraph_cartesian_product(ds::Vector{Digraph})
    @req !isempty(ds) "digraph_cartesian_product requires a non-empty list of digraphs"
    gap_list = GapObj([GapObj(d) for d in ds], recursive=true)
    return Digraph(DigraphWrap.DigraphCartesianProduct(gap_list))
end

function digraph_cartesian_product(d1::Digraph, d2::Digraph, rest::Digraph...)
    return digraph_cartesian_product(Digraph[d1, d2, rest...])
end

# ----------------------------------------------------------------------------
# 3.3-33  DigraphDirectProduct
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_direct_product(d1::Digraph, d2::Digraph) -> Digraph
    digraph_direct_product(d1::Digraph, d2::Digraph, rest::Digraph...) -> Digraph
    digraph_direct_product(ds::Vector{Digraph}) -> Digraph

Return the direct product of the given digraphs: there is an edge from
`[u, u']` to `[v, v']` iff there is an edge `u -> v` in the first factor and
an edge `u' -> v'` in the second. Vertex labels of the result are pairs of
factor labels, and the projections are available via
[`digraph_direct_product_projections`](@ref).

# Examples
```jldoctest
julia> digraph_direct_product(chain_digraph(4), cycle_digraph(3))
Digraph with 12 vertices, 9 edges
```
"""
function digraph_direct_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.DigraphDirectProduct(GapObj(d1), GapObj(d2)))
end

function digraph_direct_product(ds::Vector{Digraph})
    @req !isempty(ds) "digraph_direct_product requires a non-empty list of digraphs"
    gap_list = GapObj([GapObj(d) for d in ds], recursive=true)
    return Digraph(DigraphWrap.DigraphDirectProduct(gap_list))
end

function digraph_direct_product(d1::Digraph, d2::Digraph, rest::Digraph...)
    return digraph_direct_product(Digraph[d1, d2, rest...])
end

# ----------------------------------------------------------------------------
# 3.3-34  ConormalProduct
# ----------------------------------------------------------------------------
@doc raw"""
    conormal_product(d1::Digraph, d2::Digraph) -> Digraph

Return the conormal product of the digraphs `d1` and `d2` (both without
multiple edges): there is an edge from `[a, b]` to `[c, d]` iff there is an
edge `a -> c` in `d1` or an edge `b -> d` in `d2`.

# Examples
```jldoctest
julia> d = digraph_symmetric_closure(cycle_digraph(3));

julia> conormal_product(d, d)
Digraph with 9 vertices, 72 edges
```
"""
function conormal_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.ConormalProduct(GapObj(d1), GapObj(d2)))
end

# ----------------------------------------------------------------------------
# 3.3-35  HomomorphicProduct
# ----------------------------------------------------------------------------
@doc raw"""
    homomorphic_product(d1::Digraph, d2::Digraph) -> Digraph

Return the homomorphic product of the digraphs `d1` and `d2` (both without
multiple edges): there is an edge from `[a, b]` to `[c, d]` iff `a = c` or
there is an edge `a -> c` in `d1` and no edge `b -> d` in `d2`.

# Examples
```jldoctest
julia> d = digraph_symmetric_closure(chain_digraph(3));

julia> homomorphic_product(d, d)
Digraph with 9 vertices, 47 edges
```
"""
function homomorphic_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.HomomorphicProduct(GapObj(d1), GapObj(d2)))
end

# ----------------------------------------------------------------------------
# 3.3-36  LexicographicProduct
# ----------------------------------------------------------------------------
@doc raw"""
    lexicographic_product(d1::Digraph, d2::Digraph) -> Digraph

Return the lexicographic product of the digraphs `d1` and `d2` (both without
multiple edges): there is an edge from `[a, b]` to `[c, d]` iff there is an
edge `a -> c` in `d1`, or `a = c` and there is an edge `b -> d` in `d2`.

# Examples
```jldoctest
julia> lexicographic_product(digraph_symmetric_closure(cycle_digraph(3)),
       digraph_symmetric_closure(chain_digraph(2)))
Digraph with 6 vertices, 30 edges
```
"""
function lexicographic_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.LexicographicProduct(GapObj(d1), GapObj(d2)))
end

"""
    digraph_lexicographic_product(d1::Digraph, d2::Digraph) -> Digraph

Alias for [`lexicographic_product`](@ref).
"""
function digraph_lexicographic_product(d1::Digraph, d2::Digraph)
    return lexicographic_product(d1, d2)
end

# ----------------------------------------------------------------------------
# 3.3-37  ModularProduct
# ----------------------------------------------------------------------------
@doc raw"""
    modular_product(d1::Digraph, d2::Digraph) -> Digraph

Return the modular product of the digraphs `d1` and `d2` (both without
multiple edges): there is an edge from `[a, b]` to `[c, d]` iff `a = c`
exactly when `b = d`, and `a -> c` is an edge exactly when `b -> d` is an
edge. The complete subdigraphs of the modular product are precisely the
partial isomorphisms from `d1` to `d2`.

# Examples
```jldoctest
julia> modular_product(digraph([[1], [1, 2]]), digraph([Int[], [2]]))
Digraph with 4 vertices, 4 edges
```
"""
function modular_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.ModularProduct(GapObj(d1), GapObj(d2)))
end

"""
    digraph_modular_product(d1::Digraph, d2::Digraph) -> Digraph

Alias for [`modular_product`](@ref).
"""
function digraph_modular_product(d1::Digraph, d2::Digraph)
    return modular_product(d1, d2)
end

# ----------------------------------------------------------------------------
# 3.3-38  StrongProduct
# ----------------------------------------------------------------------------
@doc raw"""
    strong_product(d1::Digraph, d2::Digraph) -> Digraph

Return the strong product of the digraphs `d1` and `d2` (both without
multiple edges): there is an edge from `[a, b]` to `[c, d]` iff `a = c` and
`b -> d` is an edge, or `b = d` and `a -> c` is an edge, or both `a -> c` and
`b -> d` are edges.

# Examples
```jldoctest
julia> strong_product(digraph_symmetric_closure(chain_digraph(3)),
       digraph_symmetric_closure(chain_digraph(4)))
Digraph with 12 vertices, 58 edges
```
"""
function strong_product(d1::Digraph, d2::Digraph)
    return Digraph(DigraphWrap.StrongProduct(GapObj(d1), GapObj(d2)))
end

"""
    digraph_strong_product(d1::Digraph, d2::Digraph) -> Digraph

Alias for [`strong_product`](@ref).
"""
function digraph_strong_product(d1::Digraph, d2::Digraph)
    return strong_product(d1, d2)
end

# ----------------------------------------------------------------------------
# 3.3-39  DigraphCartesianProductProjections
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_cartesian_product_projections(d::Digraph) -> Vector{GapObj}

If `d` is a Cartesian product digraph, return the list of projections onto
the coordinates of the product (each projection is an idempotent
endomorphism of `d`). This attribute is only set for digraphs created with
[`digraph_cartesian_product`](@ref).

# Examples
```jldoctest
julia> d = digraph_cartesian_product(chain_digraph(3), cycle_digraph(4));

julia> length(digraph_cartesian_product_projections(d))
2
```
"""
function digraph_cartesian_product_projections(d::Digraph)
    return Vector{GapObj}(DigraphWrap.DigraphCartesianProductProjections(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-40  DigraphDirectProductProjections
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_direct_product_projections(d::Digraph) -> Vector{GapObj}

If `d` is a direct product digraph, return the list of projections onto the
coordinates of the product (each projection is an idempotent endomorphism of
`d`). This attribute is only set for digraphs created with
[`digraph_direct_product`](@ref).

# Examples
```jldoctest
julia> d = digraph_direct_product(chain_digraph(3), cycle_digraph(4));

julia> length(digraph_direct_product_projections(d))
2
```
"""
function digraph_direct_product_projections(d::Digraph)
    return Vector{GapObj}(DigraphWrap.DigraphDirectProductProjections(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-41  LineDigraph / EdgeDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    line_digraph(d::Digraph) -> Digraph

Return the line digraph of `d`: each edge of `d` becomes a vertex, and there
is an edge from the vertex of an edge `[u, v]` to the vertex of an edge
`[v, w]` for every path `u -> v -> w` in `d`. The argument is never modified.

# Examples
```jldoctest
julia> line_digraph(complete_digraph(3))
Digraph with 6 vertices, 12 edges
```
"""
function line_digraph(d::Digraph)
    return Digraph(DigraphWrap.LineDigraph(GapObj(d)))
end

"""
    edge_digraph(d::Digraph) -> Digraph

Synonym for [`line_digraph`](@ref), mirroring the GAP alias `EdgeDigraph`.
"""
function edge_digraph(d::Digraph)
    return line_digraph(d)
end

# ----------------------------------------------------------------------------
# 3.3-42  LineUndirectedDigraph / EdgeUndirectedDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    line_undirected_digraph(d::Digraph) -> Digraph

Given a symmetric digraph `d`, return the symmetric digraph whose vertices
are the edges of `d` (ignoring directions and multiplicities), with an edge
between two vertices iff the corresponding edges have a vertex in common.
The argument is never modified.

# Examples
```jldoctest
julia> line_undirected_digraph(complete_digraph(3))
Digraph with 3 vertices, 6 edges
```
"""
function line_undirected_digraph(d::Digraph)
    return Digraph(DigraphWrap.LineUndirectedDigraph(GapObj(d)))
end

"""
    edge_undirected_digraph(d::Digraph) -> Digraph

Synonym for [`line_undirected_digraph`](@ref), mirroring the GAP alias
`EdgeUndirectedDigraph`.
"""
function edge_undirected_digraph(d::Digraph)
    return line_undirected_digraph(d)
end

# ----------------------------------------------------------------------------
# 3.3-43  DoubleDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    double_digraph(d::Digraph) -> Digraph

Return the double digraph of `d`: a copy of `d` together with a duplicated
vertex set, containing the original edges, their duplicates, and the edges
`[u1, v2]` and `[u2, v1]` for every edge `[u, v]` of `d`. If `d` is mutable
it is changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> double_digraph(digraph([[2], [3], [1]]))
Digraph with 6 vertices, 12 edges
```
"""
function double_digraph(d::Digraph)
    return Digraph(DigraphWrap.DoubleDigraph(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-44  BipartiteDoubleDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    bipartite_double_digraph(d::Digraph) -> Digraph

Return the bipartite double digraph of `d`: a copy of `d` together with a
duplicated vertex set, containing only the edges `[u1, v2]` and `[u2, v1]`
for every edge `[u, v]` of `d`. The result is bipartite. If `d` is mutable it
is changed in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> bipartite_double_digraph(digraph([[2], [3], [1]]))
Digraph with 6 vertices, 6 edges
```
"""
function bipartite_double_digraph(d::Digraph)
    return Digraph(DigraphWrap.BipartiteDoubleDigraph(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-45  DigraphAddAllLoops
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_add_all_loops(d::Digraph) -> Digraph

Return the digraph obtained from `d` by adding a loop `[i, i]` at every
vertex `i` which does not already have one. If `d` is mutable the loops are
added in-place; otherwise a new immutable digraph is returned.

# Examples
```jldoctest
julia> digraph_add_all_loops(empty_digraph(5))
Digraph with 5 vertices, 5 edges
```
"""
function digraph_add_all_loops(d::Digraph)
    return Digraph(DigraphWrap.DigraphAddAllLoops(GapObj(d)))
end

"""
    digraph_add_all_loops_attr(d::Digraph) -> Digraph

Same as [`digraph_add_all_loops`](@ref), but for an immutable digraph the
result is stored as an attribute of `d`.
"""
function digraph_add_all_loops_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphAddAllLoopsAttr(GapObj(d)))
end

# ----------------------------------------------------------------------------
# 3.3-46  DistanceDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    distance_digraph(d::Digraph, i::Integer) -> Digraph
    distance_digraph(d::Digraph, list::Vector{<:Integer}) -> Digraph

Return the digraph on the same vertices as `d` in which `[u, v]` is an edge
iff the distance from `u` to `v` in `d` equals `i` (or belongs to `list`).
If `d` is mutable it is changed in-place; otherwise a new immutable digraph
is returned.

# Examples
```jldoctest
julia> d = cycle_digraph(5)
Digraph with 5 vertices, 5 edges

julia> distance_digraph(d, 2)
Digraph with 5 vertices, 5 edges

julia> distance_digraph(d, [1, 2])
Digraph with 5 vertices, 10 edges
```
"""
function distance_digraph(d::Digraph, i::Integer)
    return Digraph(DigraphWrap.DistanceDigraph(GapObj(d), Int(i)))
end

function distance_digraph(d::Digraph, list::Vector{<:Integer})
    return Digraph(DigraphWrap.DistanceDigraph(GapObj(d), _gap_verts(list)))
end

# ----------------------------------------------------------------------------
# 3.3-47  DigraphClosure
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_closure(d::Digraph, k::Integer) -> Digraph

Return the `k`-closure of the symmetric loopless digraph `d` with no multiple
edges: the unique smallest symmetric loopless digraph on the vertices of `d`
which contains all edges of `d` and in which the sum of the degrees of every
two non-adjacent vertices is less than `k`.

# Examples
```jldoctest
julia> d = digraph_remove_edges(complete_digraph(6), [[1, 2], [2, 1]])
Digraph with 6 vertices, 28 edges

julia> digraph_closure(d, 6)
Digraph with 6 vertices, 30 edges
```
"""
function digraph_closure(d::Digraph, k::Integer)
    return Digraph(DigraphWrap.DigraphClosure(GapObj(d), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.3-48  DigraphMycielskian
# ----------------------------------------------------------------------------
@doc raw"""
    digraph_mycielskian(d::Digraph) -> Digraph

Return the Mycielskian of the symmetric digraph `d`: a larger symmetric
digraph constructed from `d` whose chromatic number is one larger. If `d` is
mutable it is changed in-place; otherwise a new immutable digraph is
returned.

# Examples
```jldoctest
julia> digraph_mycielskian(cycle_digraph(2))
Digraph with 5 vertices, 10 edges
```
"""
function digraph_mycielskian(d::Digraph)
    return Digraph(DigraphWrap.DigraphMycielskian(GapObj(d)))
end

"""
    digraph_mycielskian_attr(d::Digraph) -> Digraph

Same as [`digraph_mycielskian`](@ref), but for an immutable digraph the
result is stored as an attribute of `d`.
"""
function digraph_mycielskian_attr(d::Digraph)
    return Digraph(DigraphWrap.DigraphMycielskianAttr(GapObj(d)))
end

# ############################################################################
# Section 3.4  Random digraphs
# ############################################################################

# Helper mapping the Julia `mut` symbol to the corresponding GAP filter.
# Besides mutability (`:mut` / `:immut`), structural filters can be requested
# exactly as in GAP's RandomDigraph.
function _random_filt(mut::Symbol)
    if mut === :mut
        return GAP.Globals.IsMutableDigraph
    elseif mut === :immut
        return GAP.Globals.IsImmutableDigraph
    elseif mut === :connected
        return GAP.Globals.IsConnectedDigraph
    elseif mut === :symmetric
        return GAP.Globals.IsSymmetricDigraph
    elseif mut === :acyclic
        return GAP.Globals.IsAcyclicDigraph
    elseif mut === :eulerian
        return GAP.Globals.IsEulerianDigraph
    elseif mut === :hamiltonian
        return GAP.Globals.IsHamiltonianDigraph
    else
        error("unsupported random digraph filter: $mut; use :mut, :immut, :connected, :symmetric, :acyclic, :eulerian, or :hamiltonian")
    end
end

# ----------------------------------------------------------------------------
# 3.4-1  RandomDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    random_digraph(n::Integer; mut::Symbol=:immut) -> Digraph
    random_digraph(n::Integer, p::AbstractFloat; mut::Symbol=:immut) -> Digraph

Return a random digraph with `n` vertices and no multiple edges, where each
possible edge is present with probability approximately `p` (if `p` is
omitted, a random probability is used). The keyword `mut` selects the GAP
filter: `:mut` or `:immut` for mutability, and `:connected`, `:symmetric`,
`:acyclic`, `:eulerian`, or `:hamiltonian` to force the corresponding
structure. The default is `:immut`, matching the default of all other
constructors in this package.
"""
function random_digraph(n::Integer; mut::Symbol=:immut)
    return Digraph(DigraphWrap.RandomDigraph(_random_filt(mut), Int(n)))
end

function random_digraph(n::Integer, p::AbstractFloat; mut::Symbol=:immut)
    return Digraph(DigraphWrap.RandomDigraph(_random_filt(mut), Int(n), GapObj(p)))
end

# ----------------------------------------------------------------------------
# 3.4-2  RandomMultiDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    random_multi_digraph(n::Integer) -> Digraph
    random_multi_digraph(n::Integer, m::Integer) -> Digraph

Return a random multidigraph with `n` vertices. If `m` is given, the digraph
has exactly `m` edges; otherwise the number of edges is chosen uniformly at
random from `[1 .. binomial(n, 2)]`. Multiple edges and loops may occur.
"""
function random_multi_digraph(n::Integer)
    return Digraph(DigraphWrap.RandomMultiDigraph(Int(n)))
end

function random_multi_digraph(n::Integer, m::Integer)
    return Digraph(DigraphWrap.RandomMultiDigraph(Int(n), Int(m)))
end

# ----------------------------------------------------------------------------
# 3.4-3  RandomTournament
# ----------------------------------------------------------------------------
@doc raw"""
    random_tournament(n::Integer; mut::Bool=false) -> Digraph

Return a random tournament with `n` vertices.
"""
function random_tournament(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.RandomTournament(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.4-4  RandomLattice
# ----------------------------------------------------------------------------
@doc raw"""
    random_lattice(n::Integer; mut::Bool=false) -> Digraph

Return a random lattice digraph with `m` vertices, where it is guaranteed
that `n <= m <= 2n`.
"""
function random_lattice(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.RandomLattice(_filt(mut), Int(n)))
end

# ############################################################################
# Section 3.5  Standard examples
# ############################################################################

# ----------------------------------------------------------------------------
# 3.5-1  AndrasfaiGraph
# ----------------------------------------------------------------------------
@doc raw"""
    andrasfai_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th Andrásfai graph: a triangle-free circulant graph with
`3n - 1` vertices.

# Examples
```jldoctest
julia> andrasfai_graph(4)
Digraph with 11 vertices, 44 edges
```
"""
function andrasfai_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.AndrasfaiGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-2  BananaTree
# ----------------------------------------------------------------------------
@doc raw"""
    banana_tree(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the `(n, k)`-banana tree: `n` copies of a `k`-star graph joined by
a single root vertex.

# Examples
```jldoctest
julia> banana_tree(2, 4)
Digraph with 9 vertices, 16 edges
```
"""
function banana_tree(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BananaTree(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-3  BinaryTree
# ----------------------------------------------------------------------------
@doc raw"""
    binary_tree(m::Integer; mut::Bool=false) -> Digraph

Construct the binary tree of depth `m` with `2^m - 1` vertices; all edges are
directed towards the root vertex `1`.

# Examples
```jldoctest
julia> binary_tree(4)
Digraph with 15 vertices, 14 edges
```
"""
function binary_tree(m::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BinaryTree(_filt(mut), Int(m)))
end

# ----------------------------------------------------------------------------
# 3.5-4  BinomialTreeGraph
# ----------------------------------------------------------------------------
@doc raw"""
    binomial_tree_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th binomial tree graph (an undirected tree).

# Examples
```jldoctest
julia> binomial_tree_graph(4)
Digraph with 4 vertices, 6 edges
```
"""
function binomial_tree_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BinomialTreeGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-5  BishopsGraph / BishopGraph
# ----------------------------------------------------------------------------
@doc raw"""
    bishops_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph
    bishops_graph(color::AbstractString, m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the bishop's graph of an `m` by `n` chessboard as a symmetric
digraph. The optional `color` argument is one of `"dark"`, `"light"`, or
`"both"` (the default) and restricts the result to the vertices of the
corresponding square colours.

# Examples
```jldoctest
julia> bishops_graph(4, 4)
Digraph with 16 vertices, 56 edges

julia> bishops_graph("dark", 3, 5)
Digraph with 8 vertices, 24 edges
```
"""
function bishops_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BishopsGraph(_filt(mut), Int(m), Int(n)))
end

function bishops_graph(color::AbstractString, m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BishopsGraph(_filt(mut), GapObj(String(color)), Int(m), Int(n)))
end

"""
    bishop_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph
    bishop_graph(color::AbstractString, m::Integer, n::Integer; mut::Bool=false) -> Digraph

Alias for [`bishops_graph`](@ref), mirroring the GAP alias `BishopGraph`.
"""
function bishop_graph(m::Integer, n::Integer; mut::Bool=false)
    return bishops_graph(m, n; mut=mut)
end

function bishop_graph(color::AbstractString, m::Integer, n::Integer; mut::Bool=false)
    return bishops_graph(color, m, n; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.5-6  BondyGraph
# ----------------------------------------------------------------------------
@doc raw"""
    bondy_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th Bondy graph (a hypohamiltonian graph).

# Examples
```jldoctest
julia> bondy_graph(1)
Digraph with 22 vertices, 66 edges
```
"""
function bondy_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BondyGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-7  BookGraph
# ----------------------------------------------------------------------------
@doc raw"""
    book_graph(m::Integer; mut::Bool=false) -> Digraph

Construct the `m`th book graph, the Cartesian product of a complete digraph
with two vertices and an `m + 1`-star graph.

# Examples
```jldoctest
julia> book_graph(2)
Digraph with 6 vertices, 14 edges
```
"""
function book_graph(m::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BookGraph(_filt(mut), Int(m)))
end

# ----------------------------------------------------------------------------
# 3.5-8  BurntPancakeGraph
# ----------------------------------------------------------------------------
@doc raw"""
    burnt_pancake_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th burnt pancake graph, the Cayley graph of the
hyperoctahedral group with respect to the prefix reversals.

# Examples
```jldoctest
julia> burnt_pancake_graph(3)
Digraph with 48 vertices, 144 edges
```
"""
function burnt_pancake_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.BurntPancakeGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-9  PancakeGraph
# ----------------------------------------------------------------------------
@doc raw"""
    pancake_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th pancake graph, the Cayley graph of the symmetric group
with respect to the prefix reversals.

# Examples
```jldoctest
julia> pancake_graph(4)
Digraph with 24 vertices, 72 edges
```
"""
function pancake_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.PancakeGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-10  StackedBookGraph
# ----------------------------------------------------------------------------
@doc raw"""
    stacked_book_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the `(m, n)`-stacked book graph, the Cartesian product of the
symmetric closure of the chain digraph with `n` vertices and an `m + 1`-star
graph.

# Examples
```jldoctest
julia> stacked_book_graph(1, 2)
Digraph with 4 vertices, 8 edges
```
"""
function stacked_book_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.StackedBookGraph(_filt(mut), Int(m), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-11  ChainDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    chain_digraph(n::Integer; mut::Bool=false) -> Digraph

Construct the chain digraph with `n` vertices and the edges `i -> i + 1` for
`1 <= i < n`.

# Examples
```jldoctest
julia> chain_digraph(4)
Digraph with 4 vertices, 3 edges
```
"""
function chain_digraph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.ChainDigraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-12  CirculantGraph
# ----------------------------------------------------------------------------
@doc raw"""
    circulant_graph(n::Integer, par::Vector{<:Integer}; mut::Bool=false) -> Digraph

Construct the circulant graph `Ci(n, par)`: vertex `i` is adjacent to the
vertices `i ± j` for every `j` in `par`.

# Examples
```jldoctest
julia> circulant_graph(6, [2, 3])
Digraph with 6 vertices, 18 edges
```
"""
function circulant_graph(n::Integer, par::Vector{<:Integer}; mut::Bool=false)
    return Digraph(DigraphWrap.CirculantGraph(_filt(mut), Int(n), GapObj(collect(Int, par))))
end

# ----------------------------------------------------------------------------
# 3.5-13  CompleteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    complete_digraph(n::Integer; mut::Bool=false) -> Digraph

Construct the complete digraph with `n` vertices.

# Examples
```jldoctest
julia> complete_digraph(3)
Digraph with 3 vertices, 6 edges
```
"""
function complete_digraph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.CompleteDigraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-14  CompleteBipartiteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    complete_bipartite_digraph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the complete bipartite digraph with parts of sizes `m` and `n`.

# Examples
```jldoctest
julia> complete_bipartite_digraph(2, 3)
Digraph with 5 vertices, 12 edges
```
"""
function complete_bipartite_digraph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.CompleteBipartiteDigraph(_filt(mut), Int(m), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-15  CompleteMultipartiteDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    complete_multipartite_digraph(orders::Vector{<:Integer}; mut::Bool=false) -> Digraph

Construct the complete multipartite digraph whose independent sets have the
sizes given by `orders`.

# Examples
```jldoctest
julia> complete_multipartite_digraph([2, 3])
Digraph with 5 vertices, 12 edges
```
"""
function complete_multipartite_digraph(orders::Vector{<:Integer}; mut::Bool=false)
    return Digraph(DigraphWrap.CompleteMultipartiteDigraph(_filt(mut), GapObj(collect(Int, orders))))
end

# ----------------------------------------------------------------------------
# 3.5-16  CycleDigraph / DigraphCycle
# ----------------------------------------------------------------------------
@doc raw"""
    cycle_digraph(n::Integer; mut::Bool=false) -> Digraph

Construct the cycle digraph with `n` vertices and the edges
`i -> i + 1` (mod `n`).

# Examples
```jldoctest
julia> cycle_digraph(4)
Digraph with 4 vertices, 4 edges
```
"""
function cycle_digraph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.CycleDigraph(_filt(mut), Int(n)))
end

"""
    digraph_cycle(n::Integer; mut::Bool=false) -> Digraph

Synonym for [`cycle_digraph`](@ref), mirroring the GAP alias `DigraphCycle`.
"""
function digraph_cycle(n::Integer; mut::Bool=false)
    return cycle_digraph(n; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.5-17  CycleGraph
# ----------------------------------------------------------------------------
@doc raw"""
    cycle_symmetric_digraph(n::Integer; mut::Bool=false) -> Digraph

Construct the undirected cycle graph on `n` vertices (as a symmetric
digraph), mirroring the GAP operation `CycleGraph`. The name differs from the
GAP name because Oscar already provides `cycle_graph` for the undirected
`Graphs.jl` graph with the same signature.

# Examples
```jldoctest
julia> cycle_symmetric_digraph(7)
Digraph with 7 vertices, 14 edges
```
"""
function cycle_symmetric_digraph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.CycleGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-18  EmptyDigraph / NullDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    empty_digraph(n::Integer; mut::Bool=false) -> Digraph

Construct the empty (null) digraph with `n` vertices and no edges.

# Examples
```jldoctest
julia> empty_digraph(5)
Digraph with 5 vertices, 0 edges
```
"""
function empty_digraph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.EmptyDigraph(_filt(mut), Int(n)))
end

"""
    null_digraph(n::Integer; mut::Bool=false) -> Digraph

Synonym for [`empty_digraph`](@ref), mirroring the GAP alias `NullDigraph`.
"""
function null_digraph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.NullDigraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-19  GearGraph
# ----------------------------------------------------------------------------
@doc raw"""
    gear_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th gear graph: a cycle graph on `2n` vertices with one
additional central vertex adjacent to every other vertex.

# Examples
```jldoctest
julia> gear_graph(4)
Digraph with 9 vertices, 24 edges
```
"""
function gear_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.GearGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-20  HaarGraph
# ----------------------------------------------------------------------------
@doc raw"""
    haar_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the Haar graph `H(n)`.

# Examples
```jldoctest
julia> haar_graph(4)
Digraph with 6 vertices, 6 edges
```
"""
function haar_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.HaarGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-21  HalvedCubeGraph
# ----------------------------------------------------------------------------
@doc raw"""
    halved_cube_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th halved cube graph (the graph of the `n`-demihypercube).

# Examples
```jldoctest
julia> halved_cube_graph(3)
Digraph with 4 vertices, 12 edges
```
"""
function halved_cube_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.HalvedCubeGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-22  HanoiGraph
# ----------------------------------------------------------------------------
@doc raw"""
    hanoi_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th Hanoi graph, whose vertices are the states of the Tower
of Hanoi puzzle on three towers.

# Examples
```jldoctest
julia> hanoi_graph(4)
Digraph with 81 vertices, 240 edges
```
"""
function hanoi_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.HanoiGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-23  HelmGraph
# ----------------------------------------------------------------------------
@doc raw"""
    helm_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th helm graph.

# Examples
```jldoctest
julia> helm_graph(4)
Digraph with 9 vertices, 24 edges
```
"""
function helm_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.HelmGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-24  HypercubeGraph
# ----------------------------------------------------------------------------
@doc raw"""
    hypercube_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`-dimensional hypercube graph.

# Examples
```jldoctest
julia> hypercube_graph(4)
Digraph with 16 vertices, 64 edges
```
"""
function hypercube_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.HypercubeGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-25  JohnsonDigraph
# ----------------------------------------------------------------------------
@doc raw"""
    johnson_digraph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the Johnson graph `J(n, k)`: the vertices are the `k`-subsets of
`[1 .. n]`, adjacent iff their intersection has size `k - 1`.

# Examples
```jldoctest
julia> johnson_digraph(4, 2)
Digraph with 6 vertices, 24 edges
```
"""
function johnson_digraph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.JohnsonDigraph(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-26  KellerGraph
# ----------------------------------------------------------------------------
@doc raw"""
    keller_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`-dimensional Keller graph, used for testing maximum clique
algorithms.

# Examples
```jldoctest
julia> keller_graph(3)
Digraph with 64 vertices, 2176 edges
```
"""
function keller_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.KellerGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-27  KingsGraph
# ----------------------------------------------------------------------------
@doc raw"""
    kings_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the king's graph of an `m` by `n` chessboard.

# Examples
```jldoctest
julia> kings_graph(4, 4)
Digraph with 16 vertices, 84 edges
```
"""
function kings_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.KingsGraph(_filt(mut), Int(m), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-28  KneserGraph
# ----------------------------------------------------------------------------
@doc raw"""
    kneser_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the Kneser graph `KG(n, k)`.

# Examples
```jldoctest
julia> kneser_graph(5, 2)
Digraph with 10 vertices, 30 edges
```
"""
function kneser_graph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.KneserGraph(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-29  KnightsGraph
# ----------------------------------------------------------------------------
@doc raw"""
    knights_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the knight's graph of an `m` by `n` chessboard.

# Examples
```jldoctest
julia> knights_graph(4, 4)
Digraph with 16 vertices, 48 edges
```
"""
function knights_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.KnightsGraph(_filt(mut), Int(m), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-30  LindgrenSousselierGraph
# ----------------------------------------------------------------------------
@doc raw"""
    lindgren_sousselier_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th Lindgren-Sousselier graph (a hypohamiltonian graph).

# Examples
```jldoctest
julia> lindgren_sousselier_graph(4)
Digraph with 28 vertices, 90 edges
```
"""
function lindgren_sousselier_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.LindgrenSousselierGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-31  LollipopGraph
# ----------------------------------------------------------------------------
@doc raw"""
    lollipop_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the `(m, n)`-lollipop graph.

# Examples
```jldoctest
julia> lollipop_graph(4, 4)
Digraph with 8 vertices, 20 edges
```
"""
function lollipop_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.LollipopGraph(_filt(mut), Int(m), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-32  MobiusLadderGraph
# ----------------------------------------------------------------------------
@doc raw"""
    mobius_ladder_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the Möbius ladder graph with `2n` vertices.

# Examples
```jldoctest
julia> mobius_ladder_graph(4)
Digraph with 8 vertices, 24 edges
```
"""
function mobius_ladder_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.MobiusLadderGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-33  MycielskiGraph
# ----------------------------------------------------------------------------
@doc raw"""
    mycielski_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th Mycielski graph.

# Examples
```jldoctest
julia> mycielski_graph(4)
Digraph with 11 vertices, 40 edges
```
"""
function mycielski_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.MycielskiGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-34  OddGraph
# ----------------------------------------------------------------------------
@doc raw"""
    odd_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th odd graph.

# Examples
```jldoctest
julia> odd_graph(4)
Digraph with 35 vertices, 140 edges
```
"""
function odd_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.OddGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-35  PathGraph
# ----------------------------------------------------------------------------
@doc raw"""
    path_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the path graph on `n` vertices (the symmetric closure of the chain
digraph).

# Examples
```jldoctest
julia> path_graph(4)
Digraph with 4 vertices, 6 edges
```
"""
function path_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.PathGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-36  PermutationStarGraph
# ----------------------------------------------------------------------------
@doc raw"""
    permutation_star_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the `(n, k)`-permutation star graph.

# Examples
```jldoctest
julia> permutation_star_graph(4, 3)
Digraph with 24 vertices, 72 edges
```
"""
function permutation_star_graph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.PermutationStarGraph(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-37  PetersenGraph
# ----------------------------------------------------------------------------
@doc raw"""
    petersen_digraph(; mut::Bool=false) -> Digraph

Construct the Petersen graph (10 vertices, 15 undirected edges) as a
symmetric digraph, mirroring the GAP operation `PetersenGraph`. The name
differs from the GAP name because Oscar already provides `petersen_graph`
for the undirected `Graphs.jl` graph with the same signature.

# Examples
```jldoctest
julia> petersen_digraph()
Digraph with 10 vertices, 30 edges
```
"""
function petersen_digraph(; mut::Bool=false)
    return Digraph(DigraphWrap.PetersenGraph(_filt(mut)))
end

# ----------------------------------------------------------------------------
# 3.5-38  GeneralisedPetersenGraph
# ----------------------------------------------------------------------------
@doc raw"""
    generalised_petersen_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the generalised Petersen graph `GPG(n, k)`.

# Examples
```jldoctest
julia> generalised_petersen_graph(7, 2)
Digraph with 14 vertices, 42 edges
```
"""
function generalised_petersen_graph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.GeneralisedPetersenGraph(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-39  PrismGraph
# ----------------------------------------------------------------------------
@doc raw"""
    prism_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th prism graph.

# Examples
```jldoctest
julia> prism_graph(4)
Digraph with 8 vertices, 24 edges
```
"""
function prism_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.PrismGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-40  StackedPrismGraph
# ----------------------------------------------------------------------------
@doc raw"""
    stacked_prism_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the `(n, k)`-stacked prism graph.

# Examples
```jldoctest
julia> stacked_prism_graph(5, 2)
Digraph with 10 vertices, 30 edges
```
"""
function stacked_prism_graph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.StackedPrismGraph(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-41  QueensGraph / QueenGraph
# ----------------------------------------------------------------------------
@doc raw"""
    queens_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the queen's graph of an `m` by `n` chessboard.

# Examples
```jldoctest
julia> queens_graph(4, 4)
Digraph with 16 vertices, 152 edges
```
"""
function queens_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.QueensGraph(_filt(mut), Int(m), Int(n)))
end

"""
    queen_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Alias for [`queens_graph`](@ref), mirroring the GAP alias `QueenGraph`.
"""
function queen_graph(m::Integer, n::Integer; mut::Bool=false)
    return queens_graph(m, n; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.5-42  RooksGraph / RookGraph
# ----------------------------------------------------------------------------
@doc raw"""
    rooks_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the rook's graph of an `m` by `n` chessboard.

# Examples
```jldoctest
julia> rooks_graph(4, 4)
Digraph with 16 vertices, 96 edges
```
"""
function rooks_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.RooksGraph(_filt(mut), Int(m), Int(n)))
end

"""
    rook_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Alias for [`rooks_graph`](@ref), mirroring the GAP alias `RookGraph`.
"""
function rook_graph(m::Integer, n::Integer; mut::Bool=false)
    return rooks_graph(m, n; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.5-43  SquareGridGraph / GridGraph
# ----------------------------------------------------------------------------
@doc raw"""
    square_grid_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the `n` by `k` square grid graph.

# Examples
```jldoctest
julia> square_grid_graph(4, 4)
Digraph with 16 vertices, 48 edges
```
"""
function square_grid_graph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.SquareGridGraph(_filt(mut), Int(n), Int(k)))
end

"""
    grid_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Alias for [`square_grid_graph`](@ref), mirroring the GAP alias `GridGraph`.
"""
function grid_graph(n::Integer, k::Integer; mut::Bool=false)
    return square_grid_graph(n, k; mut=mut)
end

# ----------------------------------------------------------------------------
# 3.5-44  TriangularGridGraph
# ----------------------------------------------------------------------------
@doc raw"""
    triangular_grid_graph(n::Integer, k::Integer; mut::Bool=false) -> Digraph

Construct the `n` by `k` triangular grid graph.

# Examples
```jldoctest
julia> triangular_grid_graph(3, 3)
Digraph with 9 vertices, 32 edges
```
"""
function triangular_grid_graph(n::Integer, k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.TriangularGridGraph(_filt(mut), Int(n), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-45  StarGraph
# ----------------------------------------------------------------------------
@doc raw"""
    star_graph(k::Integer; mut::Bool=false) -> Digraph

Construct the star graph with `k` vertices, in which vertex `1` is adjacent
to every other vertex.

# Examples
```jldoctest
julia> star_graph(5)
Digraph with 5 vertices, 8 edges
```
"""
function star_graph(k::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.StarGraph(_filt(mut), Int(k)))
end

# ----------------------------------------------------------------------------
# 3.5-46  TadpoleGraph
# ----------------------------------------------------------------------------
@doc raw"""
    tadpole_graph(m::Integer, n::Integer; mut::Bool=false) -> Digraph

Construct the `(m, n)`-tadpole graph: a cycle on `[1 .. m]` joined to a chain
on `[m + 1 .. m + n]`.

# Examples
```jldoctest
julia> tadpole_graph(4, 4)
Digraph with 8 vertices, 16 edges
```
"""
function tadpole_graph(m::Integer, n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.TadpoleGraph(_filt(mut), Int(m), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-47  WalshHadamardGraph
# ----------------------------------------------------------------------------
@doc raw"""
    walsh_hadamard_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the Hadamard graph built from the `n`th Walsh Hadamard matrix.

# Examples
```jldoctest
julia> walsh_hadamard_graph(4)
Digraph with 32 vertices, 256 edges
```
"""
function walsh_hadamard_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.WalshHadamardGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-48  WebGraph
# ----------------------------------------------------------------------------
@doc raw"""
    web_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th web graph.

# Examples
```jldoctest
julia> web_graph(5)
Digraph with 15 vertices, 40 edges
```
"""
function web_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.WebGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-49  WheelGraph
# ----------------------------------------------------------------------------
@doc raw"""
    wheel_graph(n::Integer; mut::Bool=false) -> Digraph

Construct the `n`th wheel graph.

# Examples
```jldoctest
julia> wheel_graph(8)
Digraph with 8 vertices, 28 edges
```
"""
function wheel_graph(n::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.WheelGraph(_filt(mut), Int(n)))
end

# ----------------------------------------------------------------------------
# 3.5-50  WindmillGraph
# ----------------------------------------------------------------------------
@doc raw"""
    windmill_graph(n::Integer, m::Integer; mut::Bool=false) -> Digraph

Construct the `(n, m)`-windmill graph: `m` copies of the complete graph on
`n` vertices sharing a single common vertex.

# Examples
```jldoctest
julia> windmill_graph(4, 3)
Digraph with 10 vertices, 36 edges
```
"""
function windmill_graph(n::Integer, m::Integer; mut::Bool=false)
    return Digraph(DigraphWrap.WindmillGraph(_filt(mut), Int(n), Int(m)))
end

# ############################################################################
# Appendix: edge-weighted digraph (not part of GAP Digraphs manual Chapter 3)
# ############################################################################
@doc raw"""
    edge_weighted_digraph(D::Digraph, weights::Vector{<:AbstractVector}) -> Digraph
    edge_weighted_digraph(adj::Vector{<:AbstractVector}, weights::Vector{<:AbstractVector}; mut::Bool=false) -> Digraph

Return an edge-weighted digraph from the digraph `D` with the given edge
`weights` (the `i`-th entry lists the weights of the edges of vertex `i` in
the order of its out-neighbours).
"""
function edge_weighted_digraph(D::Digraph, weights::Vector{<:AbstractVector})
    return Digraph(DigraphWrap.EdgeWeightedDigraph(GapObj(D), GapObj(weights, recursive=true)))
end

function edge_weighted_digraph(adj::Vector{<:AbstractVector}, weights::Vector{<:AbstractVector}; mut::Bool=false)
    d = digraph(adj; mut=mut)
    return edge_weighted_digraph(d, weights)
end
