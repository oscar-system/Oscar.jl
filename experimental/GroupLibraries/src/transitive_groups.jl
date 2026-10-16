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

`identify(L, G)` returns `(d, i)` where `d` is the degree of `G`,
and requires that `G` is transitive on `1:d`.
[`find(L::GroupLibrary, filters...)`](@ref) accepts `number_of_moved_points`
in place of `degree`; the two agree for the groups in the library.

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

################################################################################
#
#  Information
#
################################################################################

# the publications that determined the transitive groups of degree `d`
function _transitive_groups_references(d::Int)
  refs = ["TransGrp"]
  d <= 15 && push!(refs, "CHM98")
  if d <= 11
    push!(refs, "BM83")
  elseif d == 12
    push!(refs, "Roy87")
  elseif d in (14, 15)
    push!(refs, "But93")
  elseif 16 <= d <= 30
    append!(refs, ["Hul96", "Hul05"])
  elseif d == 32
    push!(refs, "CH08")
  elseif 34 <= d <= 46
    push!(refs, "HR20")
  elseif d == 48
    push!(refs, "HRT22")
  end
  return refs
end

const _TRANSITIVE_GROUPS_FILTERS = Any[order, transitivity, is_primitive]

function information(::TransitiveGroupsLibrary)
  lines = [
    "TransGrp $(_gap_package_version("transgrp")): degree at most 48;",
    "  the data for the degrees 32 and 48 are not installed with OSCAR",
  ]
  references = [
    "TransGrp", "CHM98", "BM83", "Roy87", "But93", "Hul96", "Hul05", "CH08", "HR20", "HRT22"
  ]
  return LibraryInformation("Library of transitive groups"; lines, references)
end

function information(L::TransitiveGroupsLibrary, d::IntegerUnion)
  d = Int(d)
  title = "Transitive groups of degree $d"
  has_groups(L, d) || return LibraryInformation(title; lines=["not available"])
  lines = [
    _count_line(number_of_transitive_groups(d), "group"),
    _filters_line(_TRANSITIVE_GROUPS_FILTERS),
    "identify: " * _availability_line(has_identification(L, d)),
  ]
  return LibraryInformation(
    title; lines, references=_transitive_groups_references(d), filters=_TRANSITIVE_GROUPS_FILTERS
  )
end

function information(L::TransitiveGroupsLibrary, id::Tuple{IntegerUnion,IntegerUnion})
  d, i = Int(id[1]), Int(id[2])
  G = GapObj(L[d, i])
  lines = ["group $i of the $(number_of_transitive_groups(d)) groups of degree $d"]
  GAP.Globals.HasName(G)::Bool && push!(lines, "name: " * String(GAP.Globals.Name(G)))
  return LibraryInformation(
    "Transitive group ($d, $i)"; lines, references=_transitive_groups_references(d)
  )
end
