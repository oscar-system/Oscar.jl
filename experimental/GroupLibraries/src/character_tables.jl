struct CharacterTableLibrary <: GroupLibrary{GAPGroupCharacterTable} end

@doc raw"""
    character_table_library()

Return a handle for the library of character tables,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The library consists of ordinary character tables instead of groups.
The identifier of a table is its name, a string as returned by
[`all_character_table_names`](@ref), and `L[name]` is the table
[`character_table(name)`](@ref character_table(id::String, p::Int = 0)).
As there, `name` can be any admissible name of the table;
`haskey(L, name)` tells whether `name` is admissible.
The `p`-modular Brauer table for `name` is `mod(L[name], p)`.

The library has no primary key, thus
`find(L, filters...; ordered_by = nothing)` takes any number of filters,
also none.
It supports the functions `exponent`, `is_abelian`, `is_almost_simple`,
`is_atlas_character_table`, `is_cyclic`, `is_duplicate_table`,
`is_nilpotent`, `is_perfect`, `is_quasisimple`, `is_simple`,
`is_sporadic_simple`, `is_solvable`, `is_supersolvable`,
`number_of_conjugacy_classes`, `order`.
A selection lists its tables in the order of `all_character_table_names`,
and by increasing value of `ordered_by` if that is one of these functions.

`length`, `keys` and `isempty` of a selection construct no table.
The library stores the values of the supported functions except for
`exponent`, `is_cyclic`, `is_nilpotent`, `is_solvable`, `is_supersolvable`;
with one of these as a filter, each evaluation of the selection reads all
tables of the library, which takes minutes.

[`identify`](@ref) returns the name of a table that was obtained from `L`,
and throws an exception for any other table.
`has_groups`, `has_number_of_groups` and `has_identification` are not
available.

# Examples
```jldoctest
julia> L = character_table_library()
Library of character tables

julia> S = find(L, is_sporadic_simple, !is_duplicate_table; ordered_by = order)
Selection of character tables: is_sporadic_simple, !is_duplicate_table; ordered_by = order

julia> length(S)
26

julia> println(first(S))
character table of M11

julia> collect(keys(S))[1:5]
5-element Vector{String}:
 "M11"
 "M12"
 "J1"
 "M22"
 "J2"

julia> [order(t) for t in find(L, number_of_conjugacy_classes => 5, !is_duplicate_table)]
7-element Vector{ZZRingElem}:
 20
 21
 60
 14
  8
  8
 24

julia> identify(L, L["L2(4)"])
"A5"

julia> haskey(L, "J5")
false
```
"""
character_table_library() = CharacterTableLibrary()

_name(::CharacterTableLibrary) = "character tables"
_primary_key(::CharacterTableLibrary) = nothing
_filter_attrs(::CharacterTableLibrary) = _ctbl_filter_attrs

function find(L::CharacterTableLibrary, filters...; ordered_by=nothing)
  options = (;)
  if ordered_by !== nothing
    find_index_function(ordered_by, _filter_attrs(L))
    options = (; ordered_by)
  end

  # `_translate` only rejects unsupported filters here
  return GroupLibrarySelection(L, Nothing[], filters, _translate(L, filters), options)
end

function Base.getindex(::CharacterTableLibrary, name::String)
  tbl = character_table(name)
  @req tbl !== nothing "there is no character table with name $name"
  return tbl
end

Base.haskey(::CharacterTableLibrary, name::String) = is_character_table_name(name)

# decided by identity: `character_table` caches the tables of the library
function identify(L::CharacterTableLibrary, tbl::GAPGroupCharacterTable)
  name = identifier(tbl)
  @req haskey(L, name) && L[name] === tbl "only tables from the library can be identified"
  return name
end

# `all_character_table_names` evaluates the filters without constructing
# a table in OSCAR
function _ids(S::GroupLibrarySelection{CharacterTableLibrary})
  return all_character_table_names(S.rest...; S.options...)
end

function _pairs(S::GroupLibrarySelection{CharacterTableLibrary})
  return ((name, S.library[name]) for name in _ids(S))
end

_count(S::GroupLibrarySelection{CharacterTableLibrary}) = length(_ids(S))
_isempty(S::GroupLibrarySelection{CharacterTableLibrary}) = isempty(_ids(S))

function _first(S::GroupLibrarySelection{CharacterTableLibrary})
  names = _ids(S)
  return isempty(names) ? nothing : S.library[first(names)]
end
