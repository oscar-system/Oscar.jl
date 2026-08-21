# ===========================================================================
# Chapter 7: Homomorphisms
#
# Wrappers for the operations and attributes of Chapter 7 ("Homomorphisms")
# of the GAP Digraphs manual:
#
#   7.1 Acting on digraphs
#     - 7.1-1 OnDigraphs
#     - 7.1-2 ^
#     - 7.1-3 OnMultiDigraphs
#     - 7.1-4 OnTuplesDigraphs / OnSetsDigraphs
#   7.2 Isomorphisms and canonical labellings
#     - 7.2-1 DigraphsUseNauty / DigraphsUseBliss
#     - 7.2-2 AutomorphismGroup
#     - 7.2-3 BlissAutomorphismGroup
#     - 7.2-4 NautyAutomorphismGroup
#     - 7.2-5 AutomorphismGroup (vertex colours)
#     - 7.2-6 AutomorphismGroup (vertex and edge colours)
#     - 7.2-7 BlissCanonicalLabelling / NautyCanonicalLabelling
#     - 7.2-8 BlissCanonicalLabelling / NautyCanonicalLabelling (colours)
#     - 7.2-9 BlissCanonicalDigraph / NautyCanonicalDigraph
#     - 7.2-10 DigraphGroup
#     - 7.2-11 DigraphOrbits
#     - 7.2-12 DigraphOrbitReps
#     - 7.2-13 DigraphSchreierVector
#     - 7.2-14 DigraphStabilizer
#     - 7.2-15 IsIsomorphicDigraph
#     - 7.2-16 IsIsomorphicDigraph (colours)
#     - 7.2-17 IsomorphismDigraphs
#     - 7.2-18 IsomorphismDigraphs (colours)
#     - 7.2-19 RepresentativeOutNeighbours
#     - 7.2-20 IsDigraphIsomorphism / IsDigraphAutomorphism
#     - 7.2-21 IsDigraphColouring
#     - 7.2-22 MaximalCommonSubdigraph
#     - 7.2-23 MinimalCommonSuperdigraph
#   7.3 Graph homomorphisms
#     - 7.3-1 HomomorphismDigraphsFinder
#     - 7.3-2 DigraphHomomorphism
#     - 7.3-3 HomomorphismsDigraphs / HomomorphismsDigraphsRepresentatives
#     - 7.3-4 DigraphMonomorphism
#     - 7.3-5 MonomorphismsDigraphs / MonomorphismsDigraphsRepresentatives
#     - 7.3-6 DigraphEpimorphism
#     - 7.3-7 EpimorphismsDigraphs / EpimorphismsDigraphsRepresentatives
#     - 7.3-8 DigraphEmbedding
#     - 7.3-9 EmbeddingsDigraphs / EmbeddingsDigraphsRepresentatives
#     - 7.3-10 IsDigraphHomomorphism / IsDigraphEpimorphism /
#            IsDigraphMonomorphism / IsDigraphEndomorphism
#     - 7.3-11 IsDigraphEmbedding
#     - 7.3-12 SubdigraphsMonomorphisms / SubdigraphsMonomorphismsRepresentatives
#     - 7.3-13 DigraphsRespectsColouring
#     - 7.3-14 GeneratorsOfEndomorphismMonoid / GeneratorsOfEndomorphismMonoidAttr
#     - 7.3-15 DigraphColouring
#     - 7.3-16 DigraphGreedyColouring
#     - 7.3-17 DigraphWelshPowellOrder
#     - 7.3-18 ChromaticNumber
#     - 7.3-19 DigraphCore
#     - 7.3-20 LatticeDigraphEmbedding
#     - 7.3-21 IsLatticeHomomorphism / IsLatticeEpimorphism /
#            IsLatticeEmbedding / IsLatticeMonomorphism / IsLatticeEndomorphism
#
# A permutation can be passed either as a `PermGroupElem` or as a raw GAP
# permutation (a `GapObj`); a transformation is always a GAP transformation
# (a `GapObj`). Colouring arguments (vertex or edge colours) can be given as
# Julia lists or as GAP objects; `GAP.Globals.fail` is accepted wherever the
# GAP manual allows it. Whenever a GAP operation returns `fail`, the
# corresponding Julia function raises an `ArgumentError` explaining why.
# ===========================================================================

# Convert a colouring or a list of GAP objects (e.g. a pair of permutations
# or an ordering of vertices) to the GAP object expected by the operations
# below. `GAP.Globals.fail` is itself a `GapObj` and is passed through.
_gap_colours(x) = x isa GapObj ? x : GapObj(x, recursive = true)

# Convert a permutation given as a `PermGroupElem` to the corresponding GAP
# permutation; raw GAP permutations and transformations are passed through.
_gap_perm(x::GapObj) = x
_gap_perm(x::PermGroupElem) = GapObj(x)

const _GAP_HOOK_WRAPPER = Ref{Any}(nothing)

# Convert a Julia hook `(user_param, t) -> ...` into a GAP function with two
# arguments, so that it can be passed to `HomomorphismDigraphsFinder`, whose
# GAP code checks the number of arguments of the hook.
function _gap_hook(hook::Function)
  if _GAP_HOOK_WRAPPER[] === nothing
    _GAP_HOOK_WRAPPER[] = GAP.evalstr("f -> function(user_param, t) return CallFuncList(f, [user_param, t]); end")
  end
  return _GAP_HOOK_WRAPPER[](hook)
end

# ---------------------------------------------------------------------------
# 7.1-1  OnDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    on_digraphs(D::Digraph, p::PermGroupElem) -> Digraph
    on_digraphs(D::Digraph, t::GapObj) -> Digraph

If the second argument `p` is a permutation of the vertices of `D`, return
the digraph constructed by relabelling the vertices of `D` according to `p`.
If the second argument `t` is a transformation of the vertices of `D`, return
the digraph constructed by transforming the source and range of every edge
according to `t`; a vertex not in the image of `t` becomes isolated and the
result may contain multiple edges even if `D` does not.

If `D` is mutable, the relabelling is performed directly on `D`; otherwise an
immutable copy with the vertices relabelled is returned. The vertex labels of
`D` are not retained. The same operation is available with the `^` operator.

# Examples
```jldoctest
julia> D = digraph([[3], [1, 3, 5], [1], [1, 2, 4], [2, 3, 5]]);

julia> out_neighbours(on_digraphs(D, GAP.evalstr("(1, 2)")))
5-element Vector{Vector{Int64}}:
 [2, 3, 5]
 [3]
 [2]
 [2, 1, 4]
 [1, 3, 5]

julia> out_neighbours(on_digraphs(digraph([[2], [], [2]]), GAP.Globals.Transformation(GapObj([1, 2, 1]))))
3-element Vector{Vector{Int64}}:
 [2, 2]
 []
 []
```
"""
function on_digraphs(D::Digraph, p::PermGroupElem)
  return Digraph(DigraphWrap.OnDigraphs(GapObj(D), _gap_perm(p)))
end

function on_digraphs(D::Digraph, t::GapObj)
  return Digraph(DigraphWrap.OnDigraphs(GapObj(D), t))
end

# ---------------------------------------------------------------------------
# 7.1-2  ^
# ---------------------------------------------------------------------------
@doc raw"""
    ^(D::Digraph, p) -> Digraph

Return `on_digraphs(D, p)`, where `p` is a permutation (a `PermGroupElem` or
a GAP permutation) or a transformation (a `GapObj`). This gives a short way
to apply vertex relabelling or transformations to a digraph.

# Examples
```jldoctest
julia> D = cycle_digraph(5);

julia> D ^ GAP.evalstr("(1, 5)(2, 4)") == digraph_reverse(D)
true

julia> D ^ GAP.evalstr("()") == D
true
```
"""
function Base.:^(D::Digraph, p::PermGroupElem)
  return on_digraphs(D, p)
end

function Base.:^(D::Digraph, t::GapObj)
  return on_digraphs(D, t)
end

# ---------------------------------------------------------------------------
# 7.1-3  OnMultiDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    on_multi_digraphs(D::Digraph, pair) -> Digraph
    on_multi_digraphs(D::Digraph, p1::PermGroupElem, p2::PermGroupElem) -> Digraph
    on_multi_digraphs(D::Digraph, p1::GapObj, p2::GapObj) -> Digraph

Relabel the vertices and the edges of the digraph `D` according to the
permutations `pair[1]` and `pair[2]` (respectively `p1` and `p2`). Note that
`on_digraphs(D, p) == on_multi_digraphs(D, [p, ()])`. If `D` is mutable, the
relabelling is performed directly on `D`; otherwise an immutable copy with
the vertices and edges relabelled is returned.

# Examples
```jldoctest
julia> D = cycle_digraph(3);

julia> p1 = GAP.evalstr("(1, 2)");

julia> p2 = GAP.evalstr("()");

julia> on_multi_digraphs(D, [p1, p2]) == on_digraphs(D, p1)
true

julia> on_multi_digraphs(D, p1, p2) == on_digraphs(D, p1)
true
```
"""
function on_multi_digraphs(D::Digraph, pair)
  return Digraph(DigraphWrap.OnMultiDigraphs(GapObj(D), _gap_colours(pair)))
end

function on_multi_digraphs(D::Digraph, p1::PermGroupElem, p2::PermGroupElem)
  return Digraph(DigraphWrap.OnMultiDigraphs(GapObj(D), _gap_perm(p1), _gap_perm(p2)))
end

function on_multi_digraphs(D::Digraph, p1::GapObj, p2::GapObj)
  return Digraph(DigraphWrap.OnMultiDigraphs(GapObj(D), p1, p2))
end

# ---------------------------------------------------------------------------
# 7.1-4  OnTuplesDigraphs / OnSetsDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    on_tuples_digraphs(list::Vector{Digraph}, p) -> Vector{Digraph}
    on_sets_digraphs(list::Vector{Digraph}, p) -> Vector{Digraph}

Apply the permutation `p` to a copy (with the same mutability) of each digraph
in `list` via `on_digraphs`. `on_sets_digraphs` returns the sorted output of
`on_tuples_digraphs`; if `list` is a GAP set, the result is again a set.

# Examples
```jldoctest
julia> list = [cycle_digraph(6), digraph_reverse(cycle_digraph(6))];

julia> p = GAP.evalstr("(1, 6)(2, 5)(3, 4)");

julia> on_tuples_digraphs(list, p) == reverse(list)
true

julia> on_sets_digraphs(list, p) == list
true
```
"""
function on_tuples_digraphs(list::Vector{Digraph}, p::PermGroupElem)
  result = DigraphWrap.OnTuplesDigraphs(_gap_colours([GapObj(D) for D in list]), _gap_perm(p))
  return [Digraph(result[i]) for i in 1:length(result)]
end

function on_tuples_digraphs(list::Vector{Digraph}, p::GapObj)
  result = DigraphWrap.OnTuplesDigraphs(_gap_colours([GapObj(D) for D in list]), p)
  return [Digraph(result[i]) for i in 1:length(result)]
end

@doc raw"""
    on_sets_digraphs(list::Vector{Digraph}, p) -> Vector{Digraph}

Return the sorted output of `on_tuples_digraphs(list, p)`; if `list` is a
GAP set, the result is again a set. See [`on_tuples_digraphs`](@ref).
"""
function on_sets_digraphs(list::Vector{Digraph}, p::PermGroupElem)
  result = DigraphWrap.OnSetsDigraphs(_gap_colours([GapObj(D) for D in list]), _gap_perm(p))
  return [Digraph(result[i]) for i in 1:length(result)]
end

function on_sets_digraphs(list::Vector{Digraph}, p::GapObj)
  result = DigraphWrap.OnSetsDigraphs(_gap_colours([GapObj(D) for D in list]), p)
  return [Digraph(result[i]) for i in 1:length(result)]
end

# ---------------------------------------------------------------------------
# 7.2-1  DigraphsUseNauty / DigraphsUseBliss
# ---------------------------------------------------------------------------
@doc raw"""
    digraphs_use_nauty() -> Nothing
    digraphs_use_bliss() -> Nothing

Specify whether nauty or bliss should be used by default for subsequent
computations of canonical labellings and automorphism groups. If
`NautyTracesInterface` is not available, `digraphs_use_nauty` does nothing.
These functions can be called at any point in a session; existing digraphs
remain valid, and comparisons of existing and newly created digraphs via
`is_isomorphic` and `isomorphism_digraphs` remain valid.

# Examples
```jldoctest
julia> digraphs_use_bliss()

julia> digraphs_use_nauty()
```
"""
function digraphs_use_nauty()
  return DigraphWrap.DigraphsUseNauty()
end

@doc raw"""
    digraphs_use_bliss() -> Nothing

Specify that bliss should be used by default for subsequent computations of
canonical labellings and automorphism groups; see
[`digraphs_use_nauty`](@ref).
"""
function digraphs_use_bliss()
  return DigraphWrap.DigraphsUseBliss()
end

# ---------------------------------------------------------------------------
# 7.2-2 / 7.2-5 / 7.2-6  AutomorphismGroup
# ---------------------------------------------------------------------------
@doc raw"""
    automorphism_group(D::Digraph) -> PermGroup
    automorphism_group(D::Digraph, vert_colours) -> PermGroup
    automorphism_group(D::Digraph, vert_colours, edge_colours) -> PermGroup

Return the group of automorphisms of the digraph `D`. If `D` is not a
multidigraph, the returned group is a permutation group on the vertices of
`D`; if `D` is a multidigraph, it is the direct product of a permutation
group on the vertices and a permutation group on the edges.

The colouring `vert_colours` is either a list of `nv(D)` integers (colour
`i` of vertex `i`) or a list of non-empty disjoint lists whose union is the
set of vertices of `D`; the colouring `edge_colours` is a list of `nv(D)`
lists of the same shape as `out_neighbours(D)`. Giving `GAP.Globals.fail`
for `vert_colours` [resp. `edge_colours`] is equivalent to giving all
vertices [resp. edges] the same colour. By default the group is computed
using bliss; if `NautyTracesInterface` is available it can also be computed
using nauty after calling `digraphs_use_nauty`.

# Examples
```jldoctest
julia> order(automorphism_group(complete_digraph(4)))
24

julia> order(automorphism_group(cycle_digraph(9), [[1, 4, 7], [2, 5, 8], [3, 6, 9]]))
3
```
"""
function automorphism_group(D::Digraph)
  return PermGroup(DigraphWrap.AutomorphismGroup(GapObj(D)))
end

function automorphism_group(D::Digraph, vert_colours)
  return PermGroup(DigraphWrap.AutomorphismGroup(GapObj(D), _gap_colours(vert_colours)))
end

function automorphism_group(D::Digraph, vert_colours, edge_colours)
  return PermGroup(DigraphWrap.AutomorphismGroup(GapObj(D), _gap_colours(vert_colours), _gap_colours(edge_colours)))
end

# ---------------------------------------------------------------------------
# 7.2-3  BlissAutomorphismGroup
# ---------------------------------------------------------------------------
@doc raw"""
    bliss_automorphism_group(D::Digraph) -> PermGroup
    bliss_automorphism_group(D::Digraph, vert_colours) -> PermGroup
    bliss_automorphism_group(D::Digraph, vert_colours, edge_colours) -> PermGroup

Return the automorphism group of `D` as computed by bliss. The colouring
arguments have the same form as in [`automorphism_group`](@ref); note that
for a multidigraph, an edge-coloured automorphism group can only be computed
when no two edges share the same source, range, and colour.

# Examples
```jldoctest
julia> order(bliss_automorphism_group(complete_digraph(3)))
6
```
"""
function bliss_automorphism_group(D::Digraph)
  return PermGroup(DigraphWrap.BlissAutomorphismGroup(GapObj(D)))
end

function bliss_automorphism_group(D::Digraph, vert_colours)
  return PermGroup(DigraphWrap.BlissAutomorphismGroup(GapObj(D), _gap_colours(vert_colours)))
end

function bliss_automorphism_group(D::Digraph, vert_colours, edge_colours)
  return PermGroup(DigraphWrap.BlissAutomorphismGroup(GapObj(D), _gap_colours(vert_colours), _gap_colours(edge_colours)))
end

# ---------------------------------------------------------------------------
# 7.2-4  NautyAutomorphismGroup
# ---------------------------------------------------------------------------
@doc raw"""
    nauty_automorphism_group(D::Digraph) -> PermGroup
    nauty_automorphism_group(D::Digraph, vert_colours) -> PermGroup

Return the automorphism group of `D` as computed by nauty via
`NautyTracesInterface`, which must be available (throw an `ArgumentError`
otherwise). The colouring argument has the same form as in
[`automorphism_group`](@ref).
"""
function nauty_automorphism_group(D::Digraph)
  result = DigraphWrap.NautyAutomorphismGroup(GapObj(D))
  @req result != GAP.Globals.fail "nauty_automorphism_group: the GAP package NautyTracesInterface is not available"
  return PermGroup(result)
end

function nauty_automorphism_group(D::Digraph, vert_colours)
  result = DigraphWrap.NautyAutomorphismGroup(GapObj(D), _gap_colours(vert_colours))
  @req result != GAP.Globals.fail "nauty_automorphism_group: the GAP package NautyTracesInterface is not available"
  return PermGroup(result)
end

# ---------------------------------------------------------------------------
# 7.2-7 / 7.2-8  BlissCanonicalLabelling / NautyCanonicalLabelling
# ---------------------------------------------------------------------------
@doc raw"""
    bliss_canonical_labelling(D::Digraph) -> GapObj
    bliss_canonical_labelling(D::Digraph, colours) -> GapObj
    nauty_canonical_labelling(D::Digraph) -> GapObj
    nauty_canonical_labelling(D::Digraph, colours) -> GapObj

Return the canonical labelling of `D` computed by bliss (respectively nauty).
If `D` has no multiple edges, the canonical labelling is a single permutation
of the vertices; for a multidigraph it is a pair `[p, q]` consisting of a
permutation `p` of the vertices and a permutation `q` of the edges. The
canonical representative of `D` is `on_digraphs(D, p)`, see
[`bliss_canonical_digraph`](@ref). The argument `colours` restricts the
canonical labelling to colour-preserving isomorphisms and has the same form
as the vertex colouring in [`automorphism_group`](@ref).

# Examples
```jldoctest
julia> lab = bliss_canonical_labelling(cycle_digraph(4));

julia> on_digraphs(cycle_digraph(4), lab) == bliss_canonical_digraph(cycle_digraph(4))
true
```
"""
function bliss_canonical_labelling(D::Digraph)
  return DigraphWrap.BlissCanonicalLabelling(GapObj(D))
end

function bliss_canonical_labelling(D::Digraph, colours)
  return DigraphWrap.BlissCanonicalLabelling(GapObj(D), _gap_colours(colours))
end

@doc raw"""
    nauty_canonical_labelling(D::Digraph) -> GapObj
    nauty_canonical_labelling(D::Digraph, colours) -> GapObj

Return the canonical labelling of `D` as computed by nauty; see
[`bliss_canonical_labelling`](@ref) for the meaning of the return value.
Throw an `ArgumentError` if the GAP package `NautyTracesInterface` is not
available.
"""
function nauty_canonical_labelling(D::Digraph)
  result = DigraphWrap.NautyCanonicalLabelling(GapObj(D))
  @req result != GAP.Globals.fail "nauty_canonical_labelling: the GAP package NautyTracesInterface is not available"
  return result
end

function nauty_canonical_labelling(D::Digraph, colours)
  result = DigraphWrap.NautyCanonicalLabelling(GapObj(D), _gap_colours(colours))
  @req result != GAP.Globals.fail "nauty_canonical_labelling: the GAP package NautyTracesInterface is not available"
  return result
end

# ---------------------------------------------------------------------------
# 7.2-9  BlissCanonicalDigraph / NautyCanonicalDigraph
# ---------------------------------------------------------------------------
@doc raw"""
    bliss_canonical_digraph(D::Digraph) -> Digraph
    nauty_canonical_digraph(D::Digraph) -> Digraph

Return the canonical representative of the isomorphism class of `D`, computed
by bliss (respectively nauty), i.e. `on_digraphs(D, p)` where `p` is the
canonical labelling of `D`. In particular, two digraphs are isomorphic if
and only if their canonical representatives are equal.

# Examples
```jldoctest
julia> bliss_canonical_digraph(cycle_digraph(4)) isa Digraph
true
```
"""
function bliss_canonical_digraph(D::Digraph)
  return Digraph(DigraphWrap.BlissCanonicalDigraph(GapObj(D)))
end

@doc raw"""
    nauty_canonical_digraph(D::Digraph) -> Digraph

Return the canonical representative of the isomorphism class of `D` as
computed by nauty; see [`bliss_canonical_digraph`](@ref). Throw an
`ArgumentError` if the GAP package `NautyTracesInterface` is not available.
"""
function nauty_canonical_digraph(D::Digraph)
  result = DigraphWrap.NautyCanonicalDigraph(GapObj(D))
  @req result != GAP.Globals.fail "nauty_canonical_digraph: the GAP package NautyTracesInterface is not available"
  return Digraph(result)
end

# ---------------------------------------------------------------------------
# 7.2-10  DigraphGroup
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_group(D::Digraph) -> PermGroup

Return the group of the Cayley digraph `D` if applicable. If `D` is immutable
and was created knowing a subgroup of its automorphism group, that group is
stored in the attribute `DigraphGroup` and is returned; otherwise (and always
for a mutable digraph) the full automorphism group of `D` is returned.

# Examples
```jldoctest
julia> d = cayley_digraph(symmetric_group(3));

julia> order(digraph_group(d))
6
```
"""
function digraph_group(D::Digraph)
  return PermGroup(DigraphWrap.DigraphGroup(GapObj(D)))
end

# ---------------------------------------------------------------------------
# 7.2-11 / 7.2-12 / 7.2-13 / 7.2-14  Orbits of the digraph group
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_orbits(D::Digraph) -> Vector{Vector{Int}}

Return the orbits of the action of `digraph_group(D)` on the set of vertices
of `D`.

# Examples
```jldoctest
julia> D = cayley_digraph(alternating_group(4));

julia> length(digraph_orbits(D))
1
```
"""
function digraph_orbits(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.DigraphOrbits(GapObj(D)))
end

@doc raw"""
    digraph_orbit_reps(D::Digraph) -> Vector{Int}

Return a list of orbit representatives of the action of `digraph_group(D)`
on the set of vertices of `D`.

# Examples
```jldoctest
julia> digraph_orbit_reps(cayley_digraph(alternating_group(4)))
1-element Vector{Int64}:
 1
```
"""
function digraph_orbit_reps(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphOrbitReps(GapObj(D)))
end

@doc raw"""
    digraph_schreier_vector(D::Digraph) -> Vector{Int}

Return the Schreier vector of the action of `digraph_group(D)` on the set of
vertices of `D`: a list `sch` of integers such that `sch[i] < 0` means that
`i` is an orbit representative, and `sch[i] > 0` means that
`i / gens[sch[i]]` is one step closer to the root of the Schreier tree.

# Examples
```jldoctest
julia> length(digraph_schreier_vector(cayley_digraph(alternating_group(4))))
12
```
"""
function digraph_schreier_vector(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphSchreierVector(GapObj(D)))
end

@doc raw"""
    digraph_stabilizer(D::Digraph, v::Int) -> PermGroup

Return the stabilizer of the vertex `v` under the action of
`digraph_group(D)` on the set of vertices of `D`.

# Examples
```jldoctest
julia> order(digraph_stabilizer(cayley_digraph(alternating_group(4)), 1))
1
```
"""
function digraph_stabilizer(D::Digraph, v::Int)
  return PermGroup(DigraphWrap.DigraphStabilizer(GapObj(D), v))
end

# ---------------------------------------------------------------------------
# 7.2-15 / 7.2-16  IsIsomorphicDigraph
# ---------------------------------------------------------------------------
@doc raw"""
    is_isomorphic(D1::Digraph, D2::Digraph) -> Bool
    is_isomorphic(D1::Digraph, D2::Digraph, colours1, colours2) -> Bool

Return `true` if the digraphs `D1` and `D2` are isomorphic, i.e. if there
exists a permutation `p` of the vertices of `D1` such that
`on_digraphs(D1, p) == D2`. If the colourings `colours1` and `colours2` are
given, return `true` if there exists an isomorphism that also preserves the
colourings.

# Examples
```jldoctest
julia> is_isomorphic(cycle_digraph(4), digraph_from_edges([[1, 2], [2, 3], [3, 4], [4, 1]]))
true

julia> is_isomorphic(cycle_digraph(4), cycle_digraph(4), [1, 2, 3, 4], [1, 2, 3, 4])
true
```
"""
function is_isomorphic(D1::Digraph, D2::Digraph)
  return DigraphWrap.IsIsomorphicDigraph(GapObj(D1), GapObj(D2))::Bool
end

function is_isomorphic(D1::Digraph, D2::Digraph, colours1, colours2)
  return DigraphWrap.IsIsomorphicDigraph(GapObj(D1), GapObj(D2), _gap_colours(colours1), _gap_colours(colours2))::Bool
end

# ---------------------------------------------------------------------------
# 7.2-17 / 7.2-18  IsomorphismDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    isomorphism_digraphs(D1::Digraph, D2::Digraph) -> GapObj
    isomorphism_digraphs(D1::Digraph, D2::Digraph, colours1, colours2) -> GapObj

Return a permutation `p` of the vertices of `D1` such that
`on_digraphs(D1, p) == D2`; if the colourings `colours1` and `colours2` are
given, `p` is an isomorphism of coloured digraphs. Throw an `ArgumentError`
if `D1` and `D2` are not isomorphic.

# Examples
```jldoctest
julia> D = cycle_digraph(4);

julia> p = isomorphism_digraphs(D, D);

julia> on_digraphs(D, p) == D
true
```
"""
function isomorphism_digraphs(D1::Digraph, D2::Digraph)
  result = DigraphWrap.IsomorphismDigraphs(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "isomorphism_digraphs: the digraphs are not isomorphic"
  return result
end

function isomorphism_digraphs(D1::Digraph, D2::Digraph, colours1, colours2)
  result = DigraphWrap.IsomorphismDigraphs(GapObj(D1), GapObj(D2), _gap_colours(colours1), _gap_colours(colours2))
  @req result != GAP.Globals.fail "isomorphism_digraphs: the coloured digraphs are not isomorphic"
  return result
end

# ---------------------------------------------------------------------------
# 7.2-19  RepresentativeOutNeighbours
# ---------------------------------------------------------------------------
@doc raw"""
    representative_out_neighbours(D::Digraph) -> Vector{Vector{Int}}

Return the out-neighbours of each representative of the orbits of the action
of `digraph_group(D)` on the set of vertices of `D`. If `digraph_group(D)`
is trivial, this is exactly `out_neighbours(D)`.

# Examples
```jldoctest
julia> D = digraph([[2], [3], []]);

julia> representative_out_neighbours(D) == out_neighbours(D)
true
```
"""
function representative_out_neighbours(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.RepresentativeOutNeighbours(GapObj(D)))
end

# ---------------------------------------------------------------------------
# 7.2-20  IsDigraphIsomorphism / IsDigraphAutomorphism
# ---------------------------------------------------------------------------
@doc raw"""
    is_digraph_isomorphism(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_isomorphism(src::Digraph, ran::Digraph, x, col1, col2) -> Bool
    is_digraph_automorphism(D::Digraph, x) -> Bool
    is_digraph_automorphism(D::Digraph, x, col) -> Bool

Return `true` if the permutation or transformation `x` is an isomorphism from
`src` to `ran` (respectively an automorphism of `D`), i.e. a homomorphism of
digraphs whose inverse is again a homomorphism. If colourings are given, `x`
must additionally preserve them; the colourings have the same form as in
[`automorphism_group`](@ref).

# Examples
```jldoctest
julia> D = cycle_digraph(4);

julia> p = isomorphism_digraphs(D, D);

julia> is_digraph_isomorphism(D, D, p)
true

julia> is_digraph_automorphism(D, p)
true
```
"""
function is_digraph_isomorphism(src::Digraph, ran::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphIsomorphism(GapObj(src), GapObj(ran), _gap_perm(x))::Bool
end

function is_digraph_isomorphism(src::Digraph, ran::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphIsomorphism(GapObj(src), GapObj(ran), x)::Bool
end

function is_digraph_isomorphism(src::Digraph, ran::Digraph, x, col1, col2)
  return DigraphWrap.IsDigraphIsomorphism(GapObj(src), GapObj(ran), _gap_perm(x), _gap_colours(col1), _gap_colours(col2))::Bool
end

@doc raw"""
    is_digraph_automorphism(D::Digraph, x) -> Bool
    is_digraph_automorphism(D::Digraph, x, col) -> Bool

Return `true` if the permutation or transformation `x` is an automorphism of
`D`, i.e. an isomorphism from `D` to itself; if the colouring `col` is
given, `x` must additionally preserve it. See also
[`is_digraph_isomorphism`](@ref).
"""
function is_digraph_automorphism(D::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphAutomorphism(GapObj(D), _gap_perm(x))::Bool
end

function is_digraph_automorphism(D::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphAutomorphism(GapObj(D), x)::Bool
end

function is_digraph_automorphism(D::Digraph, x, col)
  return DigraphWrap.IsDigraphAutomorphism(GapObj(D), _gap_perm(x), _gap_colours(col))::Bool
end

# ---------------------------------------------------------------------------
# 7.2-21  IsDigraphColouring
# ---------------------------------------------------------------------------
@doc raw"""
    is_digraph_colouring(D::Digraph, list::Vector{<:Integer}) -> Bool
    is_digraph_colouring(D::Digraph, t::GapObj) -> Bool

Return `true` if the list `list` (whose `i`-th entry is the colour of vertex
`i` of `D`) or the transformation `t` (whose image is the colour of each
vertex) is a proper colouring of `D`, i.e. no edge joins two vertices of the
same colour.

# Examples
```jldoctest
julia> is_digraph_colouring(cycle_digraph(4), [1, 2, 1, 2])
true
```
"""
function is_digraph_colouring(D::Digraph, list::Vector{<:Integer})
  return DigraphWrap.IsDigraphColouring(GapObj(D), _gap_colours(list))::Bool
end

function is_digraph_colouring(D::Digraph, t::GapObj)
  return DigraphWrap.IsDigraphColouring(GapObj(D), t)::Bool
end

# ---------------------------------------------------------------------------
# 7.2-22 / 7.2-23  MaximalCommonSubdigraph / MinimalCommonSuperdigraph
# ---------------------------------------------------------------------------
@doc raw"""
    maximal_common_subdigraph(D1::Digraph, D2::Digraph) -> Tuple{Digraph, GapObj, GapObj}

Return a maximal common subdigraph of the digraphs `D1` and `D2` as a triple
consisting of the subdigraph and the two transformations embedding it into
`D1` and `D2`, respectively. Throw an `ArgumentError` if `D1` and `D2` have
no common subdigraph.

# Examples
```jldoctest
julia> result = maximal_common_subdigraph(cycle_digraph(4), cycle_digraph(4));

julia> result[1] isa Digraph && length(result) == 3
true
```
"""
function maximal_common_subdigraph(D1::Digraph, D2::Digraph)
  result = DigraphWrap.MaximalCommonSubdigraph(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "maximal_common_subdigraph: the digraphs have no common subdigraph"
  return (Digraph(result[1]), result[2], result[3])
end

@doc raw"""
    minimal_common_superdigraph(D1::Digraph, D2::Digraph) -> Tuple{Digraph, GapObj, GapObj}

Return a minimal common superdigraph of the digraphs `D1` and `D2` as a
triple consisting of the superdigraph and the two transformations embedding
`D1` and `D2` into it, respectively. Throw an `ArgumentError` if `D1` and
`D2` have no common superdigraph.

# Examples
```jldoctest
julia> result = minimal_common_superdigraph(cycle_digraph(4), cycle_digraph(4));

julia> result[1] isa Digraph && length(result) == 3
true
```
"""
function minimal_common_superdigraph(D1::Digraph, D2::Digraph)
  result = DigraphWrap.MinimalCommonSuperdigraph(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "minimal_common_superdigraph: the digraphs have no common superdigraph"
  return (Digraph(result[1]), result[2], result[3])
end

# ---------------------------------------------------------------------------
# 7.3-1  HomomorphismDigraphsFinder
# ---------------------------------------------------------------------------
@doc raw"""
    homomorphism_digraphs_finder(D1::Digraph, D2::Digraph, hook, user_param,
        max_results, hint, injective, image, partial_map, colors1, colors2;
        order=nothing, aut_grp=nothing)

Find homomorphisms from the digraph `D1` to the digraph `D2` subject to the
conditions imposed by the remaining arguments, and return the argument
`user_param`.

- `hook` is a function of two arguments `(user_param, t)`, which is called
  for every homomorphism `t` found; if it returns `true`, the search stops.
  It can also be `GAP.Globals.fail`, in which case every homomorphism found
  is added to `user_param`, which must then be a mutable list.
- `max_results` is a positive integer or `GAP.Globals.infinity`; the search
  stops after `max_results` homomorphisms have been found.
- `hint` is a positive integer restricting the rank of the homomorphisms
  found, or `GAP.Globals.fail` for no restriction.
- `injective` is `0`, `1`, or `2` (or `false`/`true`): `2` means that only
  embeddings are found, `1` that only injective homomorphisms are found, and
  `0` imposes no restriction.
- `image` is a subset of the vertices of `D2`.
- `partial_map` is a partial map from `D1` to `D2` which the homomorphisms
  found must extend, or `GAP.Globals.fail`.
- `colors1` and `colors2` are colourings of `D1` and `D2` (or `fail`).
- The optional keyword arguments are `order`, an ordering of the vertices of
  `D1` for the search, and `aut_grp`, a subgroup of the automorphism group of
  `D2` up to which the homomorphisms found are unique representatives. If
  `aut_grp` is given but `order` is not, `order` defaults to `1:nv(D1)`.

# Examples
```jldoctest
julia> D = cycle_digraph(4);

julia> hook = GAP.evalstr("function(user_param, t) return false; end");

julia> homomorphism_digraphs_finder(D, D, hook, GAP.Globals.fail, 1, GAP.Globals.fail,
       false, [1, 2, 3, 4], GAP.Globals.fail, GAP.Globals.fail, GAP.Globals.fail) == GAP.Globals.fail
true
```
"""
function homomorphism_digraphs_finder(D1::Digraph, D2::Digraph, hook, user_param,
                                      max_results, hint, injective, image,
                                      partial_map, colors1, colors2;
                                      order=nothing, aut_grp=nothing)
  hook = hook isa Function ? _gap_hook(hook) : _gap_colours(hook)
  user_param = _gap_colours(user_param)
  max_results = max_results isa Integer ? max_results : _gap_colours(max_results)
  hint = hint isa Integer ? hint : _gap_colours(hint)
  injective = injective isa Bool ? Int(injective) : injective
  image = _gap_colours(image)
  partial_map = _gap_colours(partial_map)
  colors1 = _gap_colours(colors1)
  colors2 = _gap_colours(colors2)
  if aut_grp === nothing
    if order === nothing
      return DigraphWrap.HomomorphismDigraphsFinder(GapObj(D1), GapObj(D2), hook,
        user_param, max_results, hint, injective, image, partial_map, colors1,
        colors2)
    else
      return DigraphWrap.HomomorphismDigraphsFinder(GapObj(D1), GapObj(D2), hook,
        user_param, max_results, hint, injective, image, partial_map, colors1,
        colors2, _gap_colours(order))
    end
  else
    order = order === nothing ? collect(1:nv(D1)) : _gap_colours(order)
    return DigraphWrap.HomomorphismDigraphsFinder(GapObj(D1), GapObj(D2), hook,
      user_param, max_results, hint, injective, image, partial_map, colors1,
      colors2, order, _gap_colours(aut_grp))
  end
end

# ---------------------------------------------------------------------------
# 7.3-2 / 7.3-3  DigraphHomomorphism / HomomorphismsDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_homomorphism(D1::Digraph, D2::Digraph) -> GapObj

Return a homomorphism from the digraph `D1` to the digraph `D2` as a GAP
transformation, if one exists. Throw an `ArgumentError` if no homomorphism
exists.

# Examples
```jldoctest
julia> D1 = chain_digraph(3);

julia> D2 = complete_digraph(3);

julia> t = digraph_homomorphism(D1, D2);

julia> is_digraph_homomorphism(D1, D2, t)
true
```
"""
function digraph_homomorphism(D1::Digraph, D2::Digraph)
  result = DigraphWrap.DigraphHomomorphism(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "digraph_homomorphism: no homomorphism from the first digraph to the second exists"
  return result
end

@doc raw"""
    homomorphisms_digraphs(D1::Digraph, D2::Digraph) -> Vector{GapObj}
    homomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return all homomorphisms from the digraph `D1` to the digraph `D2` as GAP
transformations. `homomorphisms_digraphs_representatives` returns the
homomorphisms up to right multiplication by an element of the automorphism
group of `D2`.

# Examples
```jldoctest
julia> length(homomorphisms_digraphs(chain_digraph(2), complete_digraph(2)))
2

julia> length(homomorphisms_digraphs_representatives(chain_digraph(2), complete_digraph(2)))
1
```
"""
function homomorphisms_digraphs(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.HomomorphismsDigraphs(GapObj(D1), GapObj(D2)))
end

@doc raw"""
    homomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return the homomorphisms from `D1` to `D2` up to right multiplication by an
element of the automorphism group of `D2`; see
[`homomorphisms_digraphs`](@ref).
"""
function homomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.HomomorphismsDigraphsRepresentatives(GapObj(D1), GapObj(D2)))
end

# ---------------------------------------------------------------------------
# 7.3-4 / 7.3-5  DigraphMonomorphism / MonomorphismsDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_monomorphism(D1::Digraph, D2::Digraph) -> GapObj

Return an injective homomorphism (a monomorphism) from the digraph `D1` to
the digraph `D2` as a GAP transformation, if one exists. Throw an
`ArgumentError` if no monomorphism exists.

# Examples
```jldoctest
julia> D1 = chain_digraph(2);

julia> D2 = complete_digraph(2);

julia> t = digraph_monomorphism(D1, D2);

julia> is_digraph_monomorphism(D1, D2, t)
true
```
"""
function digraph_monomorphism(D1::Digraph, D2::Digraph)
  result = DigraphWrap.DigraphMonomorphism(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "digraph_monomorphism: no monomorphism from the first digraph to the second exists"
  return result
end

@doc raw"""
    monomorphisms_digraphs(D1::Digraph, D2::Digraph) -> Vector{GapObj}
    monomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return all injective homomorphisms (monomorphisms) from the digraph `D1` to
the digraph `D2` as GAP transformations; the `_representatives` variant
returns them up to right multiplication by an element of the automorphism
group of `D2`.

# Examples
```jldoctest
julia> length(monomorphisms_digraphs(chain_digraph(2), complete_digraph(2)))
2
```
"""
function monomorphisms_digraphs(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.MonomorphismsDigraphs(GapObj(D1), GapObj(D2)))
end

@doc raw"""
    monomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return the monomorphisms from `D1` to `D2` up to right multiplication by an
element of the automorphism group of `D2`; see
[`monomorphisms_digraphs`](@ref).
"""
function monomorphisms_digraphs_representatives(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.MonomorphismsDigraphsRepresentatives(GapObj(D1), GapObj(D2)))
end

# ---------------------------------------------------------------------------
# 7.3-6 / 7.3-7  DigraphEpimorphism / EpimorphismsDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_epimorphism(D1::Digraph, D2::Digraph) -> GapObj

Return a surjective homomorphism (an epimorphism) from the digraph `D1` to
the digraph `D2` as a GAP transformation, if one exists. Throw an
`ArgumentError` if no epimorphism exists.

# Examples
```jldoctest
julia> D = complete_digraph(2);

julia> t = digraph_epimorphism(D, D);

julia> is_digraph_epimorphism(D, D, t)
true
```
"""
function digraph_epimorphism(D1::Digraph, D2::Digraph)
  result = DigraphWrap.DigraphEpimorphism(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "digraph_epimorphism: no epimorphism from the first digraph to the second exists"
  return result
end

@doc raw"""
    epimorphisms_digraphs(D1::Digraph, D2::Digraph) -> Vector{GapObj}
    epimorphisms_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return all surjective homomorphisms (epimorphisms) from the digraph `D1` to
the digraph `D2` as GAP transformations; the `_representatives` variant
returns them up to right multiplication by an element of the automorphism
group of `D2`.

# Examples
```jldoctest
julia> length(epimorphisms_digraphs(complete_digraph(2), complete_digraph(2)))
2
```
"""
function epimorphisms_digraphs(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.EpimorphismsDigraphs(GapObj(D1), GapObj(D2)))
end

@doc raw"""
    epimorphisms_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return the epimorphisms from `D1` to `D2` up to right multiplication by an
element of the automorphism group of `D2`; see
[`epimorphisms_digraphs`](@ref).
"""
function epimorphisms_digraphs_representatives(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.EpimorphismsDigraphsRepresentatives(GapObj(D1), GapObj(D2)))
end

# ---------------------------------------------------------------------------
# 7.3-8 / 7.3-9  DigraphEmbedding / EmbeddingsDigraphs
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_embedding(D1::Digraph, D2::Digraph) -> GapObj

Return an embedding of the digraph `D1` into the digraph `D2` as a GAP
transformation, if one exists. An embedding is a monomorphism which maps
non-edges of `D1` to non-edges of `D2`, i.e. an isomorphism from `D1` to the
subdigraph of `D2` induced by the image. Throw an `ArgumentError` if no
embedding exists.

# Examples
```jldoctest
julia> D1 = chain_digraph(2);

julia> D2 = chain_digraph(3);

julia> t = digraph_embedding(D1, D2);

julia> is_digraph_embedding(D1, D2, t)
true
```
"""
function digraph_embedding(D1::Digraph, D2::Digraph)
  result = DigraphWrap.DigraphEmbedding(GapObj(D1), GapObj(D2))
  @req result != GAP.Globals.fail "digraph_embedding: no embedding of the first digraph into the second exists"
  return result
end

@doc raw"""
    embeddings_digraphs(D1::Digraph, D2::Digraph) -> Vector{GapObj}
    embeddings_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return all embeddings of the digraph `D1` into the digraph `D2` as GAP
transformations; the `_representatives` variant returns them up to right
multiplication by an element of the automorphism group of `D2`.

# Examples
```jldoctest
julia> length(embeddings_digraphs(chain_digraph(2), chain_digraph(3)))
2
```
"""
function embeddings_digraphs(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.EmbeddingsDigraphs(GapObj(D1), GapObj(D2)))
end

@doc raw"""
    embeddings_digraphs_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return the embeddings of `D1` into `D2` up to right multiplication by an
element of the automorphism group of `D2`; see
[`embeddings_digraphs`](@ref).
"""
function embeddings_digraphs_representatives(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.EmbeddingsDigraphsRepresentatives(GapObj(D1), GapObj(D2)))
end

# ---------------------------------------------------------------------------
# 7.3-10  IsDigraphHomomorphism / IsDigraphEpimorphism / IsDigraphMonomorphism /
#         IsDigraphEndomorphism
# ---------------------------------------------------------------------------
@doc raw"""
    is_digraph_homomorphism(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_homomorphism(src::Digraph, ran::Digraph, x, col1, col2) -> Bool
    is_digraph_epimorphism(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_epimorphism(src::Digraph, ran::Digraph, x, col1, col2) -> Bool
    is_digraph_monomorphism(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_monomorphism(src::Digraph, ran::Digraph, x, col1, col2) -> Bool
    is_digraph_endomorphism(D::Digraph, x) -> Bool
    is_digraph_endomorphism(D::Digraph, x, col) -> Bool

Return `true` if the permutation or transformation `x` is a homomorphism,
epimorphism, monomorphism, or endomorphism, respectively, of the digraphs
involved; if colourings are given, `x` must additionally preserve them.
An endomorphism of `D` is a homomorphism from `D` to itself.

# Examples
```jldoctest
julia> D1 = chain_digraph(3);

julia> D2 = complete_digraph(3);

julia> t = digraph_homomorphism(D1, D2);

julia> is_digraph_homomorphism(D1, D2, t)
true

julia> is_digraph_endomorphism(D2, GAP.Globals.Transformation(GapObj([1, 1, 1])))
false
```
"""
function is_digraph_homomorphism(src::Digraph, ran::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphHomomorphism(GapObj(src), GapObj(ran), _gap_perm(x))::Bool
end

function is_digraph_homomorphism(src::Digraph, ran::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphHomomorphism(GapObj(src), GapObj(ran), x)::Bool
end

function is_digraph_homomorphism(src::Digraph, ran::Digraph, x, col1, col2)
  return DigraphWrap.IsDigraphHomomorphism(GapObj(src), GapObj(ran), _gap_perm(x), _gap_colours(col1), _gap_colours(col2))::Bool
end

@doc raw"""
    is_digraph_epimorphism(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_epimorphism(src::Digraph, ran::Digraph, x, col1, col2) -> Bool

Return `true` if the permutation or transformation `x` is a surjective
homomorphism (an epimorphism) from `src` to `ran`; if colourings are given,
`x` must additionally preserve them. See also
[`is_digraph_homomorphism`](@ref).
"""
function is_digraph_epimorphism(src::Digraph, ran::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphEpimorphism(GapObj(src), GapObj(ran), _gap_perm(x))::Bool
end

function is_digraph_epimorphism(src::Digraph, ran::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphEpimorphism(GapObj(src), GapObj(ran), x)::Bool
end

function is_digraph_epimorphism(src::Digraph, ran::Digraph, x, col1, col2)
  return DigraphWrap.IsDigraphEpimorphism(GapObj(src), GapObj(ran), _gap_perm(x), _gap_colours(col1), _gap_colours(col2))::Bool
end

@doc raw"""
    is_digraph_monomorphism(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_monomorphism(src::Digraph, ran::Digraph, x, col1, col2) -> Bool

Return `true` if the permutation or transformation `x` is an injective
homomorphism (a monomorphism) from `src` to `ran`; if colourings are given,
`x` must additionally preserve them. See also
[`is_digraph_homomorphism`](@ref).
"""
function is_digraph_monomorphism(src::Digraph, ran::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphMonomorphism(GapObj(src), GapObj(ran), _gap_perm(x))::Bool
end

function is_digraph_monomorphism(src::Digraph, ran::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphMonomorphism(GapObj(src), GapObj(ran), x)::Bool
end

function is_digraph_monomorphism(src::Digraph, ran::Digraph, x, col1, col2)
  return DigraphWrap.IsDigraphMonomorphism(GapObj(src), GapObj(ran), _gap_perm(x), _gap_colours(col1), _gap_colours(col2))::Bool
end

@doc raw"""
    is_digraph_endomorphism(D::Digraph, x) -> Bool
    is_digraph_endomorphism(D::Digraph, x, col) -> Bool

Return `true` if the permutation or transformation `x` is an endomorphism of
`D`, i.e. a homomorphism from `D` to itself; if the colouring `col` is
given, `x` must additionally preserve it. See also
[`is_digraph_homomorphism`](@ref).
"""
function is_digraph_endomorphism(D::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphEndomorphism(GapObj(D), _gap_perm(x))::Bool
end

function is_digraph_endomorphism(D::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphEndomorphism(GapObj(D), x)::Bool
end

function is_digraph_endomorphism(D::Digraph, x, col)
  return DigraphWrap.IsDigraphEndomorphism(GapObj(D), _gap_perm(x), _gap_colours(col))::Bool
end

# ---------------------------------------------------------------------------
# 7.3-11  IsDigraphEmbedding
# ---------------------------------------------------------------------------
@doc raw"""
    is_digraph_embedding(src::Digraph, ran::Digraph, x) -> Bool
    is_digraph_embedding(src::Digraph, ran::Digraph, x, col1, col2) -> Bool

Return `true` if the permutation or transformation `x` is an embedding of
`src` into `ran`, i.e. a monomorphism such that the inverse of `x` is a
monomorphism from the subdigraph of `ran` induced by the image of `x` to
`src`. If colourings are given, `x` must additionally preserve them.

# Examples
```jldoctest
julia> D1 = chain_digraph(2);

julia> D2 = chain_digraph(3);

julia> t = digraph_embedding(D1, D2);

julia> is_digraph_embedding(D1, D2, t)
true
```
"""
function is_digraph_embedding(src::Digraph, ran::Digraph, x::PermGroupElem)
  return DigraphWrap.IsDigraphEmbedding(GapObj(src), GapObj(ran), _gap_perm(x))::Bool
end

function is_digraph_embedding(src::Digraph, ran::Digraph, x::GapObj)
  return DigraphWrap.IsDigraphEmbedding(GapObj(src), GapObj(ran), x)::Bool
end

function is_digraph_embedding(src::Digraph, ran::Digraph, x, col1, col2)
  return DigraphWrap.IsDigraphEmbedding(GapObj(src), GapObj(ran), _gap_perm(x), _gap_colours(col1), _gap_colours(col2))::Bool
end

# ---------------------------------------------------------------------------
# 7.3-12  SubdigraphsMonomorphisms / SubdigraphsMonomorphismsRepresentatives
# ---------------------------------------------------------------------------
@doc raw"""
    subdigraphs_monomorphisms(D1::Digraph, D2::Digraph) -> Vector{GapObj}
    subdigraphs_monomorphisms_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return all injective homomorphisms from the digraph `D1` to the digraph `D2`
whose images define the (not necessarily induced) subdigraphs of `D2` that
are isomorphic to `D1`; the `_representatives` variant returns them up to
right multiplication by an element of the automorphism group of `D2`.

# Examples
```jldoctest
julia> length(subdigraphs_monomorphisms(chain_digraph(2), chain_digraph(3)))
2
```
"""
function subdigraphs_monomorphisms(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.SubdigraphsMonomorphisms(GapObj(D1), GapObj(D2)))
end

@doc raw"""
    subdigraphs_monomorphisms_representatives(D1::Digraph, D2::Digraph) -> Vector{GapObj}

Return the monomorphisms found by `subdigraphs_monomorphisms(D1, D2)` up to
right multiplication by an element of the automorphism group of `D2`; see
[`subdigraphs_monomorphisms`](@ref).
"""
function subdigraphs_monomorphisms_representatives(D1::Digraph, D2::Digraph)
  return Vector{GapObj}(DigraphWrap.SubdigraphsMonomorphismsRepresentatives(GapObj(D1), GapObj(D2)))
end

# ---------------------------------------------------------------------------
# 7.3-13  DigraphsRespectsColouring
# ---------------------------------------------------------------------------
@doc raw"""
    digraphs_respects_colouring(src::Digraph, ran::Digraph, x, col1, col2) -> Bool

Return `true` if the permutation or transformation `x` maps the colouring
`col1` of `src` to the colouring `col2` of `ran`, i.e. if
`col1[i] == col2[i ^ x]` for every vertex `i` of `src`. The colourings have
the same form as in [`automorphism_group`](@ref).

# Examples
```jldoctest
julia> D = cycle_digraph(4);

julia> id = GAP.Globals.Transformation(GapObj([1, 2, 3, 4]));

julia> digraphs_respects_colouring(D, D, id, [1, 2, 3, 4], [1, 2, 3, 4])
true
```
"""
function digraphs_respects_colouring(src::Digraph, ran::Digraph, x, col1, col2)
  return DigraphWrap.DigraphsRespectsColouring(GapObj(src), GapObj(ran), _gap_perm(x), _gap_colours(col1), _gap_colours(col2))::Bool
end

# ---------------------------------------------------------------------------
# 7.3-14  GeneratorsOfEndomorphismMonoid / GeneratorsOfEndomorphismMonoidAttr
# ---------------------------------------------------------------------------
@doc raw"""
    generators_of_endomorphism_monoid(D::Digraph) -> Vector{GapObj}
    generators_of_endomorphism_monoid(D::Digraph, colors) -> Vector{GapObj}
    generators_of_endomorphism_monoid(D::Digraph, colors, limit::Int) -> Vector{GapObj}
    generators_of_endomorphism_monoid_attr(D::Digraph) -> Vector{GapObj}

Return a generating set of the monoid of endomorphisms of the digraph `D`
as a list of GAP transformations. If `colors` is given, only
colour-preserving endomorphisms are considered; the optional argument
`limit` bounds the computation. `generators_of_endomorphism_monoid_attr`
stores the result of `generators_of_endomorphism_monoid(D)` in an attribute
of an immutable digraph, so that it is not recomputed on future calls.

# Examples
```jldoctest
julia> !isempty(generators_of_endomorphism_monoid(complete_digraph(3)))
true
```
"""
function generators_of_endomorphism_monoid(D::Digraph)
  return Vector{GapObj}(DigraphWrap.GeneratorsOfEndomorphismMonoid(GapObj(D)))
end

function generators_of_endomorphism_monoid(D::Digraph, colors)
  return Vector{GapObj}(DigraphWrap.GeneratorsOfEndomorphismMonoid(GapObj(D), _gap_colours(colors)))
end

function generators_of_endomorphism_monoid(D::Digraph, colors, limit::Int)
  return Vector{GapObj}(DigraphWrap.GeneratorsOfEndomorphismMonoid(GapObj(D), _gap_colours(colors), limit))
end

@doc raw"""
    generators_of_endomorphism_monoid_attr(D::Digraph) -> Vector{GapObj}

Return the same value as `generators_of_endomorphism_monoid(D)`, but store it
in an attribute of the immutable digraph `D`, so that it is not recomputed on
future calls; see [`generators_of_endomorphism_monoid`](@ref).
"""
function generators_of_endomorphism_monoid_attr(D::Digraph)
  return Vector{GapObj}(DigraphWrap.GeneratorsOfEndomorphismMonoidAttr(GapObj(D)))
end

# ---------------------------------------------------------------------------
# 7.3-15  DigraphColouring
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_colouring(D::Digraph, n::Int) -> GapObj

Return a proper `n`-colouring of the digraph `D` as a GAP transformation,
if one exists. A proper colouring is a labelling of the vertices such that
adjacent vertices have different labels; note that a digraph with loops has
no proper colouring. Throw an `ArgumentError` if `D` has no proper
`n`-colouring.

# Examples
```jldoctest
julia> D = cycle_digraph(4);

julia> t = digraph_colouring(D, 2);

julia> is_digraph_colouring(D, t)
true
```
"""
function digraph_colouring(D::Digraph, n::Int)
  result = DigraphWrap.DigraphColouring(GapObj(D), n)
  @req result != GAP.Globals.fail "digraph_colouring: the digraph has no proper n-colouring"
  return result
end

# ---------------------------------------------------------------------------
# 7.3-16  DigraphGreedyColouring
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_greedy_colouring(D::Digraph) -> GapObj
    digraph_greedy_colouring(D::Digraph, order::Vector{<:Integer}) -> GapObj
    digraph_greedy_colouring(D::Digraph, func::Function) -> GapObj

Return a greedy proper colouring of the digraph `D` as a GAP transformation.
If the optional second argument is a list `order`, the vertices of `D` are
coloured greedily in the order given by `order`; if it is a function `func`,
then `func` is called with the digraph `D` and must return an ordering of
its vertices. If no second argument is given, the order
`digraph_welsh_powell_order(D)` is used. Throw an `ArgumentError` if `D` has
no proper colouring, e.g. because it has loops.

# Examples
```jldoctest
julia> D = cycle_digraph(4);

julia> is_digraph_colouring(D, digraph_greedy_colouring(D))
true

julia> is_digraph_colouring(D, digraph_greedy_colouring(D, [1, 3, 2, 4]))
true
```
"""
function digraph_greedy_colouring(D::Digraph)
  result = DigraphWrap.DigraphGreedyColouring(GapObj(D))
  @req result != GAP.Globals.fail "digraph_greedy_colouring: the digraph has no proper colouring"
  return result
end

function digraph_greedy_colouring(D::Digraph, order::Vector{<:Integer})
  result = DigraphWrap.DigraphGreedyColouring(GapObj(D), _gap_colours(order))
  @req result != GAP.Globals.fail "digraph_greedy_colouring: the digraph has no proper colouring"
  return result
end

function digraph_greedy_colouring(D::Digraph, func::Function)
  return digraph_greedy_colouring(D, func(D))
end

# ---------------------------------------------------------------------------
# 7.3-17  DigraphWelshPowellOrder
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_welsh_powell_order(D::Digraph) -> Vector{Int}

Return an ordering of the vertices of `D` from highest to lowest degree,
which can be used by [`digraph_greedy_colouring`](@ref).

# Examples
```jldoctest
julia> length(digraph_welsh_powell_order(cycle_digraph(4)))
4
```
"""
function digraph_welsh_powell_order(D::Digraph)
  return Vector{Int}(DigraphWrap.DigraphWelshPowellOrder(GapObj(D)))
end

# ---------------------------------------------------------------------------
# 7.3-18  ChromaticNumber
# ---------------------------------------------------------------------------
@doc raw"""
    chromatic_number(D::Digraph) -> Int

Return the chromatic number of the digraph `D`, i.e. the least non-negative
integer `n` such that there is a proper colouring of `D` with `n` colours.

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

# ---------------------------------------------------------------------------
# 7.3-19  DigraphCore
# ---------------------------------------------------------------------------
@doc raw"""
    digraph_core(D::Digraph) -> Vector{Int}

Return a list of vertices corresponding to the core of `D`, i.e. the minimal
subdigraph of `D` which is a homomorphic image of `D`. In particular, the
subdigraph of `D` induced by this list is isomorphic to the core of `D`.

# Examples
```jldoctest
julia> digraph_core(digraph_symmetric_closure(cycle_digraph(8)))
2-element Vector{Int64}:
 1
 2
```
"""
function digraph_core(D::Digraph)
  return Vector{Int64}(DigraphWrap.DigraphCore(GapObj(D)))
end

# ---------------------------------------------------------------------------
# 7.3-20  LatticeDigraphEmbedding
# ---------------------------------------------------------------------------
@doc raw"""
    lattice_digraph_embedding(L1::Digraph, L2::Digraph) -> GapObj

If `L1` and `L2` are lattice digraphs, return an injective homomorphism from
`L1` to `L2` which is a lattice homomorphism, as a GAP transformation. Throw
an `ArgumentError` if no such homomorphism exists.

# Examples
```jldoctest
julia> L = digraph([[1, 2, 3, 4], [2, 3, 4], [3, 4], [4]]);

julia> t = lattice_digraph_embedding(L, L);

julia> is_lattice_homomorphism(L, L, t)
true
```
"""
function lattice_digraph_embedding(L1::Digraph, L2::Digraph)
  result = DigraphWrap.LatticeDigraphEmbedding(GapObj(L1), GapObj(L2))
  @req result != GAP.Globals.fail "lattice_digraph_embedding: no embedding of the first lattice digraph into the second exists"
  return result
end

# ---------------------------------------------------------------------------
# 7.3-21  IsLatticeHomomorphism / IsLatticeEpimorphism / IsLatticeEmbedding /
#         IsLatticeMonomorphism / IsLatticeEndomorphism
# ---------------------------------------------------------------------------
@doc raw"""
    is_lattice_homomorphism(L1::Digraph, L2::Digraph, map) -> Bool
    is_lattice_epimorphism(L1::Digraph, L2::Digraph, map) -> Bool
    is_lattice_embedding(L1::Digraph, L2::Digraph, map) -> Bool
    is_lattice_monomorphism(L1::Digraph, L2::Digraph, map) -> Bool
    is_lattice_endomorphism(L::Digraph, map) -> Bool

Return `true` if the transformation `map` is a lattice homomorphism,
epimorphism, embedding, monomorphism, or endomorphism, respectively, of the
lattice digraphs involved. A lattice homomorphism is a digraph homomorphism
which also preserves the meet and the join.

# Examples
```jldoctest
julia> L = digraph([[1, 2, 3, 4], [2, 3, 4], [3, 4], [4]]);

julia> t = lattice_digraph_embedding(L, L);

julia> is_lattice_homomorphism(L, L, t) && is_lattice_endomorphism(L, t)
true
```
"""
function is_lattice_homomorphism(L1::Digraph, L2::Digraph, map::GapObj)
  return DigraphWrap.IsLatticeHomomorphism(GapObj(L1), GapObj(L2), map)::Bool
end

@doc raw"""
    is_lattice_epimorphism(L1::Digraph, L2::Digraph, map) -> Bool

Return `true` if the transformation `map` is a surjective lattice
homomorphism (a lattice epimorphism) from `L1` to `L2`; see
[`is_lattice_homomorphism`](@ref).
"""
function is_lattice_epimorphism(L1::Digraph, L2::Digraph, map::GapObj)
  return DigraphWrap.IsLatticeEpimorphism(GapObj(L1), GapObj(L2), map)::Bool
end

@doc raw"""
    is_lattice_embedding(L1::Digraph, L2::Digraph, map) -> Bool

Return `true` if the transformation `map` is a lattice embedding from `L1`
to `L2`; see [`is_lattice_homomorphism`](@ref).
"""
function is_lattice_embedding(L1::Digraph, L2::Digraph, map::GapObj)
  return DigraphWrap.IsLatticeEmbedding(GapObj(L1), GapObj(L2), map)::Bool
end

@doc raw"""
    is_lattice_monomorphism(L1::Digraph, L2::Digraph, map) -> Bool

Return `true` if the transformation `map` is an injective lattice
homomorphism (a lattice monomorphism) from `L1` to `L2`; see
[`is_lattice_homomorphism`](@ref).
"""
function is_lattice_monomorphism(L1::Digraph, L2::Digraph, map::GapObj)
  return DigraphWrap.IsLatticeMonomorphism(GapObj(L1), GapObj(L2), map)::Bool
end

@doc raw"""
    is_lattice_endomorphism(L::Digraph, map) -> Bool

Return `true` if the transformation `map` is a lattice endomorphism of `L`,
i.e. a lattice homomorphism from `L` to itself; see
[`is_lattice_homomorphism`](@ref).
"""
function is_lattice_endomorphism(L::Digraph, map::GapObj)
  return DigraphWrap.IsLatticeEndomorphism(GapObj(L), map)::Bool
end
