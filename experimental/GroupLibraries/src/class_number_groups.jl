struct GroupsWithClassNumberLibrary{T} <: GroupLibrary{T} end

@doc raw"""
    groups_with_class_number_library(::Type{T} = GAPGroup)

Return a handle for the library of groups with few conjugacy classes,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The primary key is `number_of_conjugacy_classes`; the group `L[k, i]` is the
`i`-th group with `k` conjugacy classes, up to isomorphism.
The groups are provided by the GAP package
`SmallClassNr` [SmallClassNr](@cite), for `k` up to 14.

A group is returned as `PcGroup` if it is solvable and as `PermGroup`
otherwise; if `T` is given, each group is converted to type `T`.

[`find(L::GroupLibrary, filters...)`](@ref) supports the same functions as
for [`small_groups_library`](@ref).
The library knows the `order` of its groups without constructing them.

# Examples
```jldoctest
julia> L = groups_with_class_number_library()
Library of groups with given class number

julia> L[5, 6]
Pc group of order 21

julia> length(find(L, 14))
93

julia> [describe(G) for G in find(L, 6, order => 1:20)]
5-element Vector{String}:
 "C6"
 "C3 : C4"
 "D12"
 "D18"
 "(C3 x C3) : C2"
```
"""
function groups_with_class_number_library(::Type{T}=GAPGroup) where {T}
  return GroupsWithClassNumberLibrary{T}()
end

_name(::GroupsWithClassNumberLibrary) = "groups with given class number"
_primary_key(::GroupsWithClassNumberLibrary) = number_of_conjugacy_classes
_key_type(::GroupsWithClassNumberLibrary) = ZZRingElem
_filter_attrs(::GroupsWithClassNumberLibrary) = _group_filter_attrs
_number(::GroupsWithClassNumberLibrary, k) = number_of_groups_with_class_number(k)

function Base.getindex(
  ::GroupsWithClassNumberLibrary{GAPGroup}, k::IntegerUnion, i::IntegerUnion
)
  return group_with_class_number(k, i)
end

function Base.getindex(
  ::GroupsWithClassNumberLibrary{T}, k::IntegerUnion, i::IntegerUnion
) where {T}
  return group_with_class_number(T, k, i)
end

function identify(::GroupsWithClassNumberLibrary, G::GAPGroup)
  return group_with_class_number_identification(G)
end

function has_groups(::GroupsWithClassNumberLibrary, k::IntegerUnion)
  return has_groups_with_class_number(k)
end

function has_number_of_groups(::GroupsWithClassNumberLibrary, k::IntegerUnion)
  return has_number_of_groups_with_class_number(k)
end

function has_identification(::GroupsWithClassNumberLibrary, k::IntegerUnion)
  return has_groups_with_class_number_identification(k)
end

function _gap_selection(S::GroupLibrarySelection{<:GroupsWithClassNumberLibrary}, k)
  return (GAP.Globals.NrConjugacyClasses, GAP.Obj(k), S.gaprest...)
end

function _pairs(S::GroupLibrarySelection{<:GroupsWithClassNumberLibrary}, k)
  isempty(S.rest) && return _pairs_by_index(S, k)

  lib = S.library
  _check_groups(lib, k)
  it = GAP.Globals.IteratorSmallClassNrGroups(_gap_selection(S, k)...)
  return (
    (Tuple{ZZRingElem,ZZRingElem}(GAP.Globals.IdClassNr(G)), _wrap(lib, G)) for
    G in _ConsumingGapIterator(it)
  )
end

function _count(S::GroupLibrarySelection{<:GroupsWithClassNumberLibrary}, k)
  isempty(S.rest) && return _number(S.library, k)

  _check_groups(S.library, k)
  # answered from tables if the order is the only further filter
  return GAP.Globals.NrSmallClassNrGroups(_gap_selection(S, k)...)
end
