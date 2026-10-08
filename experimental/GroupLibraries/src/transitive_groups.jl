struct TransitiveGroupsLibrary <: GroupLibrary{PermGroup} end

@doc raw"""
    transitive_groups_library()

Return a handle for the library of transitive permutation groups,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The primary key is the `degree`; the group `L[d, i]` is the `i`-th
transitive group on `d` points, up to permutation isomorphism.
The groups are provided by the GAP package `TransGrp` [TransGrp](@cite),
for degrees up to 48; the data for the degrees 32 and 48 are not installed
with OSCAR. The numbering for degree up to 15 is that of [CHM98](@cite).

A permutation group is regarded as a group on its moved points:
`identify(L, G)` returns `(d, i)` where `d` is `number_of_moved_points(G)`,
which can be smaller than `degree(G)`.
The two agree for the groups in the library, and
[`find(L::GroupLibrary, filters...)`](@ref) accepts `number_of_moved_points`
in place of `degree`.

[`find(L::GroupLibrary, filters...)`](@ref) supports the functions
`degree`, `exponent`, `is_abelian`, `is_almost_simple`, `is_cyclic`,
`is_nilpotent`, `is_perfect`, `is_primitive`, `is_quasisimple`, `is_simple`,
`is_sporadic_simple`, `is_solvable`, `is_supersolvable`, `is_transitive`,
`number_of_conjugacy_classes`, `number_of_moved_points`, `order`,
`transitivity`.
A selection with filters besides the degree constructs the matching groups
of one degree together.

# Examples
```jldoctest
julia> L = transitive_groups_library()
Library of transitive groups

julia> L[5, 4]
Alternating group of degree 5

julia> collect(find(L, degree => 3:5, is_abelian))
4-element Vector{PermGroup}:
 Alternating group of degree 3
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 5

julia> length(find(L, 30))
5712
```
"""
transitive_groups_library() = TransitiveGroupsLibrary()

_name(::TransitiveGroupsLibrary) = "transitive groups"
_primary_key(::TransitiveGroupsLibrary) = degree
_key_type(::TransitiveGroupsLibrary) = Int
_filter_attrs(::TransitiveGroupsLibrary) = _permgroup_filter_attrs
_number(::TransitiveGroupsLibrary, d) = number_of_transitive_groups(Int(d))

function _is_primary_key(::TransitiveGroupsLibrary, f)
  return f === degree || f === number_of_moved_points
end

function Base.getindex(::TransitiveGroupsLibrary, d::IntegerUnion, i::IntegerUnion)
  return transitive_group(Int(d), Int(i))
end

identify(::TransitiveGroupsLibrary, G::PermGroup) = transitive_group_identification(G)

has_groups(::TransitiveGroupsLibrary, d::IntegerUnion) = has_transitive_groups(Int(d))

function has_number_of_groups(::TransitiveGroupsLibrary, d::IntegerUnion)
  return has_number_of_transitive_groups(Int(d))
end

function has_identification(::TransitiveGroupsLibrary, d::IntegerUnion)
  return has_transitive_group_identification(Int(d))
end

# GAP's library starts at degree 2 and filters only whole degrees
_uses_gap_selection(S::GroupLibrarySelection, d) = d > 1 && !isempty(S.rest)

function _pairs(S::GroupLibrarySelection{TransitiveGroupsLibrary}, d)
  _uses_gap_selection(S, d) || return _pairs_by_index(S, d)

  _check_groups(S.library, d)
  K = GAP.Globals.AllTransitiveGroups(GAP.Globals.NrMovedPoints, d, S.gaprest...)
  return (((d, GAP.Globals.TransitiveIdentification(G)::Int), PermGroup(G, d)) for G in K)
end

function _first(S::GroupLibrarySelection{TransitiveGroupsLibrary}, d)
  _uses_gap_selection(S, d) || return _first_by_iteration(S, d)

  _check_groups(S.library, d)
  G = GAP.Globals.OneTransitiveGroup(GAP.Globals.NrMovedPoints, d, S.gaprest...)
  return G === GAP.Globals.fail ? nothing : PermGroup(G, d)
end
