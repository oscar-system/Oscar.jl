struct SmallGroupsLibrary{T} <: GroupLibrary{T} end

@doc raw"""
    small_groups_library(::Type{T} = GAPGroup)

Return a handle for the library of groups of small order,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The primary key is the `order`; the group `L[n, i]` is the `i`-th group of
order `n`, up to isomorphism.
The groups are provided by the GAP packages `SmallGrp` [SmallGrp](@cite),
`SOTGrps` [SOTGrps](@cite) and `SglPPow` [SglPPow](@cite).

A group is returned as `PcGroup` if it is solvable and as `PermGroup`
otherwise; if `T` is given, each group is converted to type `T`.

[`find(L::GroupLibrary, filters...)`](@ref) supports the functions
`exponent`, `is_abelian`, `is_almost_simple`, `is_cyclic`, `is_nilpotent`,
`is_perfect`, `is_quasisimple`, `is_simple`, `is_sporadic_simple`,
`is_solvable`, `is_supersolvable`, `number_of_conjugacy_classes`, `order`.
For most orders the library knows `is_abelian`, `is_nilpotent`, `is_solvable`
and `is_supersolvable` of its groups without constructing them.

# Examples
```jldoctest
julia> L = small_groups_library()
Library of small groups

julia> L[60, 5]
Permutation group of degree 5 and order 60

julia> length(find(L, order => 512, !is_abelian))
10494183

julia> small_groups_library(PermGroup)[8, 3]
Permutation group of degree 4 and order 8
```
"""
small_groups_library(::Type{T}=GAPGroup) where {T} = SmallGroupsLibrary{T}()

_name(::SmallGroupsLibrary) = "small groups"
_primary_key(::SmallGroupsLibrary) = order
_key_type(::SmallGroupsLibrary) = ZZRingElem
_filter_attrs(::SmallGroupsLibrary) = _group_filter_attrs
_number(::SmallGroupsLibrary, n) = number_of_small_groups(n)

function Base.getindex(::SmallGroupsLibrary{GAPGroup}, n::IntegerUnion, i::IntegerUnion)
  return small_group(n, i)
end

function Base.getindex(::SmallGroupsLibrary{T}, n::IntegerUnion, i::IntegerUnion) where {T}
  return small_group(T, n, i)
end

identify(::SmallGroupsLibrary, G::GAPGroup) = small_group_identification(G)

has_groups(::SmallGroupsLibrary, n::IntegerUnion) = has_small_groups(n)
has_number_of_groups(::SmallGroupsLibrary, n::IntegerUnion) = has_number_of_small_groups(n)
has_identification(::SmallGroupsLibrary, n::IntegerUnion) = has_small_group_identification(n)

function _pairs(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  _check_groups(S.library, n)
  return ((id, S.library[id]) for id in _ids(S, n))
end

function _ids(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  isempty(S.rest) && return ((n, ZZRingElem(i)) for i in 1:Int(_number(S.library, n)))

  _check_groups(S.library, n)
  ids = GAP.Globals.IdsOfAllSmallGroups(GAP.Globals.Size, GAP.Obj(n), S.gaprest...)
  return ((n, ZZRingElem(id[2])) for id in ids)
end

function _count(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  isempty(S.rest) && return _number(S.library, n)

  _check_groups(S.library, n)
  # answered from tables where the library indexes the given properties
  return GAP.Globals.NumberSmallGroups(GAP.Globals.Size, GAP.Obj(n), S.gaprest...)
end

function _first(S::GroupLibrarySelection{<:SmallGroupsLibrary}, n)
  _check_groups(S.library, n)
  G = GAP.Globals.OneSmallGroup(GAP.Globals.Size, GAP.Obj(n), S.gaprest...)
  return G === GAP.Globals.fail ? nothing : _wrap(S.library, G)
end
