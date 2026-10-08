struct PrimitiveGroupsLibrary <: GroupLibrary{PermGroup} end

@doc raw"""
    primitive_groups_library()

Return a handle for the library of primitive permutation groups,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The primary key is the `degree`; the group `L[d, i]` is the `i`-th
primitive group on `d` points, up to permutation isomorphism.
The groups are provided by the GAP package `PrimGrp` [PrimGrp](@cite).

A permutation group is regarded as a group on its moved points:
`identify(L, G)` returns `(d, i)` where `d` is `number_of_moved_points(G)`,
which can be smaller than `degree(G)`.
The two agree for the groups in the library, and
[`find(L::GroupLibrary, filters...)`](@ref) accepts `number_of_moved_points`
in place of `degree`.

[`find(L::GroupLibrary, filters...)`](@ref) supports the same functions as
for [`transitive_groups_library`](@ref).

# Examples
```jldoctest
julia> L = primitive_groups_library()
Library of primitive groups

julia> L[10, 1]
Permutation group of degree 10 and order 60

julia> G = stabilizer(symmetric_group(5), 1)[1];

julia> degree(G), identify(L, G)
(5, (4, 2))

julia> S = find(L, degree => 10, !is_solvable, order => 1:1000)
Selection of primitive groups: degree => 10, !is_solvable, order => 1:1000

julia> collect(keys(S))
6-element Vector{Tuple{Int64, Int64}}:
 (10, 1)
 (10, 2)
 (10, 3)
 (10, 4)
 (10, 5)
 (10, 6)
```
"""
primitive_groups_library() = PrimitiveGroupsLibrary()

_name(::PrimitiveGroupsLibrary) = "primitive groups"
_primary_key(::PrimitiveGroupsLibrary) = degree
_key_type(::PrimitiveGroupsLibrary) = Int
_filter_attrs(::PrimitiveGroupsLibrary) = _permgroup_filter_attrs
_number(::PrimitiveGroupsLibrary, d) = number_of_primitive_groups(Int(d))

function _is_primary_key(::PrimitiveGroupsLibrary, f)
  return f === degree || f === number_of_moved_points
end

function Base.getindex(::PrimitiveGroupsLibrary, d::IntegerUnion, i::IntegerUnion)
  return primitive_group(Int(d), Int(i))
end

identify(::PrimitiveGroupsLibrary, G::PermGroup) = primitive_group_identification(G)

has_groups(::PrimitiveGroupsLibrary, d::IntegerUnion) = has_primitive_groups(Int(d))

function has_number_of_groups(::PrimitiveGroupsLibrary, d::IntegerUnion)
  return has_number_of_primitive_groups(Int(d))
end

function has_identification(::PrimitiveGroupsLibrary, d::IntegerUnion)
  return has_primitive_group_identification(Int(d))
end

function _pairs(S::GroupLibrarySelection{PrimitiveGroupsLibrary}, d)
  isempty(S.rest) && return _pairs_by_index(S, d)

  _check_groups(S.library, d)
  it = GAP.Globals.PrimitiveGroupsIterator(GAP.Globals.NrMovedPoints, d, S.gaprest...)
  return (
    ((d, GAP.Globals.PrimitiveIdentification(G)::Int), PermGroup(G, d)) for
    G in _ConsumingGapIterator(it)
  )
end
