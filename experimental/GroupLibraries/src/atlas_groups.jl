struct AtlasGroupsLibrary{T} <: GroupLibrary{T} end

const _AnyAtlasGroup = Union{PermGroup,MatGroup}

@doc raw"""
    atlas_groups_library()
    atlas_groups_library(::Type{T}) where T <: Union{PermGroup, MatGroup}

Return a handle for the Atlas of Group Representations [ATLAS](@cite),
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The library consists of representations of groups.
The primary key is the name of the group, a string as accepted by
[`atlas_group`](@ref); [`find(L::GroupLibrary, filters...)`](@ref) takes it
as leading argument or as `:name => name`.
The identifier of a representation is a dictionary as returned by
[`all_atlas_group_infos`](@ref), and `L[info]` is the group `atlas_group(info)`.
The representations of a group are not numbered, since the contents of the
Atlas are not fixed; hence there is no `L[name, i]`.

Without `T` the library consists of all representations.
If `T` is `PermGroup` or `MatGroup` then it consists of the permutation
representations or the matrix representations, respectively.
In each case `first(find(L, name))` is the group that [`atlas_group`](@ref)
returns for `name` and the same `T`.

[`find(L::GroupLibrary, filters...)`](@ref) supports the functions
`degree`, `is_primitive`, `is_transitive`, `rank_action`, `transitivity`
for permutation groups and
`base_ring`, `character`, `characteristic`, `dim` for matrix groups.
They are decided from the table of contents of the Atlas;
the data of a group are fetched only when the group is constructed.

`has_number_of_groups(L, name)` tells whether `name` is the name of a group in
the Atlas, and `has_groups(L, name)` whether `L` contains a representation
of it.
[`identify`](@ref) is not available.

# Examples
```jldoctest
julia> L = atlas_groups_library()
Library of Atlas groups

julia> S = find(L, "A5", degree => [5, 6])
Selection of Atlas groups: name => "A5", degree => [5, 6]

julia> collect(keys(S))
2-element Vector{Dict{Symbol, Any}}:
 Dict(:constituents => [1, 4], :repname => "A5G1-p5B0", :degree => 5, :name => "A5")
 Dict(:constituents => [1, 5], :repname => "A5G1-p6B0", :degree => 6, :name => "A5")

julia> first(S)
Permutation group of degree 5 and order 60

julia> length(find(L, "A5")), length(find(atlas_groups_library(MatGroup), "A5"))
(18, 15)

julia> info = only(keys(find(L, "A5", dim => 4, characteristic => 3)));

julia> L[info]
Matrix group of degree 4
  over prime field of characteristic 3
```
"""
atlas_groups_library() = AtlasGroupsLibrary{_AnyAtlasGroup}()
atlas_groups_library(::Type{PermGroup}) = AtlasGroupsLibrary{PermGroup}()
atlas_groups_library(::Type{MatGroup}) = AtlasGroupsLibrary{MatGroup}()

_name(::AtlasGroupsLibrary) = "Atlas groups"
_primary_key(::AtlasGroupsLibrary) = :name
_key_type(::AtlasGroupsLibrary) = String
_filter_attrs(::AtlasGroupsLibrary) = _atlas_group_filter_attrs

_number(::AtlasGroupsLibrary{_AnyAtlasGroup}, name) = number_of_atlas_groups(name)
_number(::AtlasGroupsLibrary{T}, name) where {T} = number_of_atlas_groups(T, name)

# whether the representation described by `info` belongs to the library
_contains(::AtlasGroupsLibrary{_AnyAtlasGroup}, info::Dict) = true
_contains(::AtlasGroupsLibrary{PermGroup}, info::Dict) = haskey(info, :degree)
_contains(::AtlasGroupsLibrary{MatGroup}, info::Dict) = haskey(info, :dim)

# the entry for `name` in the table of contents of the Atlas, or `fail`
_toc_entry(name::String) = GAP.Globals.AGR.InfoForName(GapObj(name))::GapObj

function _key_values(::AtlasGroupsLibrary, x)
  x isa String && return [x]
  x isa AbstractVector{String} && return x
  return nothing
end

# several names can denote the same group, e.g. "A5" and "L2(4)"
function _canonical_name(name::String)
  entry = _toc_entry(name)
  return entry === GAP.Globals.fail ? name : String(entry[1])
end

function _normalize_keys(::AtlasGroupsLibrary, names)
  return sort!(unique!(String[_canonical_name(name) for name in names]))
end

# `all_atlas_group_infos` translates the filters itself, when the selection
# is evaluated; here `find` only rejects what is not a supported function
function _translate(L::AtlasGroupsLibrary, rest::Tuple)
  for f in rest
    f isa Pair || f isa Function ||
      throw(ArgumentError("expected a function or a pair, got $f"))
    find_index_function(f isa Pair ? f[1] : f, _filter_attrs(L))
  end
  return Any[]
end

function Base.getindex(L::AtlasGroupsLibrary{T}, info::Dict) where {T}
  @req _contains(L, info) "$(info[:repname]) does not describe a group of type $T"
  return atlas_group(info)
end

function identify(L::AtlasGroupsLibrary, G::GAPGroup)
  throw(ArgumentError("identification is not available for $(_name(L))"))
end

has_groups(L::AtlasGroupsLibrary, name::String) = _number(L, name) > 0

function has_number_of_groups(::AtlasGroupsLibrary, name::String)
  return _toc_entry(name) !== GAP.Globals.fail
end

has_identification(::AtlasGroupsLibrary, name::String) = false

function _ids(S::GroupLibrarySelection{<:AtlasGroupsLibrary}, name)
  infos = all_atlas_group_infos(name, S.rest...)
  return filter!(info -> _contains(S.library, info), infos)
end

function _pairs(S::GroupLibrarySelection{<:AtlasGroupsLibrary}, name)
  return ((info, atlas_group(info)) for info in _ids(S, name))
end

function _count(S::GroupLibrarySelection{<:AtlasGroupsLibrary}, name)
  isempty(S.rest) && return _number(S.library, name)
  return length(_ids(S, name))
end

# decided without fetching a group
_isempty(S::GroupLibrarySelection{<:AtlasGroupsLibrary}, name) = _count(S, name) == 0
