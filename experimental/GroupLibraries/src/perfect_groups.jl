struct PerfectGroupsLibrary{T} <: GroupLibrary{T} end

@doc raw"""
    perfect_groups_library(::Type{T} = PermGroup)

Return a handle for the library of perfect groups,
see [`Oscar.GroupLibrary`](@ref) for what one can do with it.

The primary key is the `order`; the group `L[n, i]` is the `i`-th perfect
group of order `n`, up to isomorphism.
The library contains the perfect groups of order less than $2 \cdot 10^6$,
see [HP89](@cite) and [Hul22](@cite);
they are provided by the GAP package `PerfGrp` [PerfGrp](@cite).

The groups have type `T`, which can be `PermGroup` or `FPGroup`.

[`find(L::GroupLibrary, filters...)`](@ref) supports every function that is
defined for groups of type `T`.
The library stores no properties of its groups, hence a filter besides
the order is evaluated on each group of the given orders.

# Examples
```jldoctest
julia> L = perfect_groups_library()
Library of perfect groups

julia> L[120, 1]
Permutation group of degree 24 and order 120

julia> perfect_groups_library(FPGroup)[120, 1]
Finitely presented group of order 120

julia> collect(find(L, order => 1:200, !is_simple))
2-element Vector{PermGroup}:
 Permutation group of degree 1 and order 1
 Permutation group of degree 24 and order 120
```
"""
function perfect_groups_library(::Type{T}=PermGroup) where {T<:Union{PermGroup,FPGroup}}
  return PerfectGroupsLibrary{T}()
end

_name(::PerfectGroupsLibrary) = "perfect groups"
_primary_key(::PerfectGroupsLibrary) = order
_key_type(::PerfectGroupsLibrary) = Int
_filter_attrs(::PerfectGroupsLibrary) = nothing
_number(::PerfectGroupsLibrary, n) = number_of_perfect_groups(n)

function Base.getindex(::PerfectGroupsLibrary{T}, n::IntegerUnion, i::IntegerUnion) where {T}
  return perfect_group(T, n, i)
end

identify(::PerfectGroupsLibrary, G::GAPGroup) = perfect_group_identification(G)

has_groups(::PerfectGroupsLibrary, n::IntegerUnion) = has_perfect_groups(Int(n))

function has_number_of_groups(::PerfectGroupsLibrary, n::IntegerUnion)
  return has_number_of_perfect_groups(Int(n))
end

function has_identification(::PerfectGroupsLibrary, n::IntegerUnion)
  return has_perfect_group_identification(Int(n))
end

################################################################################
#
#  Information
#
################################################################################

# the publications that determined the perfect groups of order `n`
_perfect_groups_references(n::Int) = ["PerfGrp", n < 10^6 ? "HP89" : "Hul22"]

function information(::PerfectGroupsLibrary)
  lines = ["PerfGrp $(_gap_package_version("perfgrp")): order less than 2 * 10^6"]
  return LibraryInformation(
    "Library of perfect groups"; lines, references=["PerfGrp", "HP89", "Hul22"]
  )
end

function information(L::PerfectGroupsLibrary, n::IntegerUnion)
  n = Int(n)
  title = "Perfect groups of order $n"
  has_groups(L, n) || return LibraryInformation(title; lines=["not available"])
  lines = [
    _count_line(number_of_perfect_groups(n), "group"),
    "identify: " * _availability_line(has_identification(L, n)),
  ]
  return LibraryInformation(title; lines, references=_perfect_groups_references(n))
end

function information(L::PerfectGroupsLibrary, id::Tuple{IntegerUnion,IntegerUnion})
  n, i = Int(id[1]), Int(id[2])
  @req has_groups(L, n) "perfect groups of order $n are not available"
  N = number_of_perfect_groups(n)
  @req 1 <= i <= N "there are only $N perfect groups of order $n, not $i"
  lines = ["group $i of the $N groups of order $n"]
  return LibraryInformation(
    "Perfect group ($n, $i)"; lines, references=_perfect_groups_references(n)
  )
end
