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

################################################################################
#
#  Information
#
################################################################################

# the publications that determined the groups with `k` conjugacy classes
function _groups_with_class_number_references(k::Int)
  refs = ["SmallClassNr"]
  if k <= 5
    append!(refs, ["Mil11", "Bur11"])
  elseif k <= 7
    push!(refs, "Pol68")
  elseif k == 8
    push!(refs, "Kos74")
  elseif k == 9
    push!(refs, "OS76")
  elseif k <= 11
    push!(refs, "V-LV-L85")
  elseif k == 12
    push!(refs, "V-LV-L86")
  elseif k <= 14
    push!(refs, "VS07")
  end
  return refs
end

const _GROUPS_WITH_CLASS_NUMBER_FILTERS = Any[order]

function information(::GroupsWithClassNumberLibrary)
  lines = [
    "SmallClassNr $(_gap_package_version("smallclassnr")): at most 14 conjugacy classes",
  ]
  references = [
    "SmallClassNr", "Lan03", "Mil11", "Bur11", "Pol68", "Kos74", "OS76", "V-LV-L85",
    "V-LV-L86", "VS07",
  ]
  return LibraryInformation("Library of groups with given class number"; lines, references)
end

function information(L::GroupsWithClassNumberLibrary, k::IntegerUnion)
  k = Int(k)
  title = "Groups with $k conjugacy classes"
  has_groups(L, k) || return LibraryInformation(title; lines=["not available"])
  lines = [
    _count_line(number_of_groups_with_class_number(k), "group"),
    _filters_line(_GROUPS_WITH_CLASS_NUMBER_FILTERS),
    "identify: " * _availability_line(has_identification(L, k)),
  ]
  return LibraryInformation(
    title;
    lines,
    references=_groups_with_class_number_references(k),
    filters=_GROUPS_WITH_CLASS_NUMBER_FILTERS,
  )
end

# the identifiers of a group in the other libraries, as GAP's SmallClassNr knows them
const _CLASS_NUMBER_CONVERSIONS = [
  ("small group", :IdClassNrToIdGroup),
  ("perfect group", :IdClassNrToPerfGrp),
  ("primitive group", :IdClassNrToPrimGrp),
  ("transitive group", :IdClassNrToTransGrp),
]

function information(L::GroupsWithClassNumberLibrary, id::Tuple{IntegerUnion,IntegerUnion})
  k, i = Int(id[1]), Int(id[2])
  @req has_groups(L, k) "groups with $k conjugacy classes are not available"
  N = number_of_groups_with_class_number(k)
  @req 1 <= i <= N "there are only $N groups with $k classes, not $i"
  lines = ["group $i of the $N groups with $k conjugacy classes"]
  for (name, f) in _CLASS_NUMBER_CONVERSIONS
    x = getproperty(GAP.Globals, f)(k, i)
    x === GAP.Globals.fail || push!(lines, "$name ($(x[1]), $(x[2]))")
  end
  x = GAP.Globals.IdClassNrToAtlasName(k, i)
  x === GAP.Globals.fail || push!(lines, "Atlas group " * String(x))
  return LibraryInformation(
    "Group ($k, $i) with $k conjugacy classes";
    lines,
    references=_groups_with_class_number_references(k),
  )
end
