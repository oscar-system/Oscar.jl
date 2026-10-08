################################################################################
#
#  Libraries of groups behind one interface
#
#    L = small_groups_library()      handle, a subtype of `GroupLibrary`
#    L[n, i]                         one group
#    S = find(L, filters...)         lazy selection, a `GroupLibrarySelection`
#    iterate(S), length(S), keys(S)  groups, their number, their ids
#
#  A selection is evaluated one value of the primary key (order, degree, ...)
#  at a time. A library implements the functions listed under "Protocol".
#
################################################################################

@doc raw"""
    GroupLibrary{T}

Abstract supertype of the handles for libraries of groups; `T` is the type of
the groups a handle returns.

A handle `L` is created by one of
[`small_groups_library`](@ref),
[`transitive_groups_library`](@ref),
[`primitive_groups_library`](@ref),
[`perfect_groups_library`](@ref),
[`groups_with_class_number_library`](@ref).
Each group in a library has an identifier. For the libraries above it is a
pair `(k, i)` where `k` is the value of the *primary key* of the library
(the order, the degree, ...) and `i` numbers the groups with that value.

- `L[k, i]` and `L[(k, i)]` return the group with identifier `(k, i)`.
- [`find(L::GroupLibrary, filters...)`](@ref) selects the groups with given
  properties.
- [`identify`](@ref) returns the identifier of a given group.
- [`has_groups`](@ref), [`has_number_of_groups`](@ref) and
  [`has_identification`](@ref) tell which part of a library is available.

# Examples
```jldoctest
julia> L = transitive_groups_library()
Library of transitive groups

julia> G = L[4, 3]
Permutation group of degree 4

julia> identify(L, G)
(4, 3)
```
"""
abstract type GroupLibrary{T} end

# The groups of `library` whose primary key lies in `keys` and which satisfy
# the filters `rest`; `gaprest` is `rest` in the form GAP's selection
# functions expect.
struct GroupLibrarySelection{L<:GroupLibrary,K<:IntegerUnion}
  library::L
  keys::Vector{K}
  rest::Tuple
  gaprest::Vector{Any}
end

################################################################################
#
#  Protocol
#
################################################################################

# Required:
#   _name(L)::String           used for printing
#   _primary_key(L)            the function whose value is the primary key
#   _key_type(L)               type of the two entries of an identifier
#   _filter_attrs(L)           supported filters, or `nothing` for any function
#   _number(L, k)              number of groups with primary key `k`
#   getindex, identify, has_groups, has_number_of_groups, has_identification
#
# Optional, for a selection `S` and a value `k` of the primary key:
#   _pairs(S, k)               lazy iterator of `(id, group)`
#   _ids(S, k)                 lazy iterator of ids
#   _count(S, k)               number of groups
#   _first(S, k)               first group, or `nothing`

_is_primary_key(L::GroupLibrary, f) = f === _primary_key(L)

_wrap(::GroupLibrary{GAPGroup}, G::GapObj) = _oscar_group(G)
_wrap(::GroupLibrary{T}, G::GapObj) where {T} = T(_oscar_group(G))

function _translate(L::GroupLibrary, rest::Tuple)
  attrs = _filter_attrs(L)
  attrs === nothing || return translate_group_library_args(rest; filter_attrs=attrs)

  for f in rest
    f isa Pair || f isa Function ||
      throw(ArgumentError("expected a function or a pair, got $f"))
  end
  return Any[]
end

function _check_groups(L::GroupLibrary, k)
  @req has_groups(L, k) "$(_name(L)) with $(_primary_key(L)) $k are not available"
end

# all groups with primary key `k` are constructed and tested one by one
function _pairs_by_index(S::GroupLibrarySelection{L,K}, k) where {L,K}
  lib = S.library
  _check_groups(lib, k)
  candidates = (((k, K(i)), lib[k, i]) for i in 1:Int(_number(lib, k)))
  isempty(S.rest) && return candidates

  return Iterators.filter(p -> _matches_group_library_filters(p[2], S.rest), candidates)
end

_pairs(S::GroupLibrarySelection, k) = _pairs_by_index(S, k)

function _ids(S::GroupLibrarySelection{L,K}, k) where {L,K}
  isempty(S.rest) || return (id for (id, _) in _pairs(S, k))

  # no group is needed to list all identifiers
  return ((k, K(i)) for i in 1:Int(_number(S.library, k)))
end

function _count(S::GroupLibrarySelection, k)
  isempty(S.rest) && return _number(S.library, k)
  return count(Returns(true), _ids(S, k))
end

function _first_by_iteration(S::GroupLibrarySelection, k)
  next = iterate(_pairs(S, k))
  return next === nothing ? nothing : next[1][2]
end

_first(S::GroupLibrarySelection, k) = _first_by_iteration(S, k)

# A GAP iterator as a Julia iterator. Each one here is created to be run
# through once, so it is consumed rather than copied first, as GAP.jl's own
# iteration does; `PrimitiveGroupsIterator` cannot be copied in PrimGrp 4.0.3.
struct _ConsumingGapIterator
  it::GapObj
end

Base.IteratorSize(::Type{_ConsumingGapIterator}) = Base.SizeUnknown()

function Base.iterate(x::_ConsumingGapIterator, _=nothing)
  GAPWrap.IsDoneIterator(x.it) && return nothing
  return GAPWrap.NextIterator(x.it), nothing
end

################################################################################
#
#  User interface
#
################################################################################

@doc raw"""
    find(L::GroupLibrary, filters...)

Return a lazy selection of the groups in `L` that satisfy all of `filters`.
Each filter has one of the following forms.

- `func => value` selects groups for which `func` returns `value`
- `func => list` selects groups for which `func` returns an element of `list`
- `func` selects groups for which `func` returns `true`
- `!func` selects groups for which `func` returns `false`

At least one filter must restrict the primary key of `L`.
A leading integer or vector of integers abbreviates such a filter,
so `find(L, 16)` means `find(L, order => 16)` for a library whose primary
key is the order.
The documentation of each library lists the functions it supports.

Creating the selection `S` computes nothing. Afterwards

- iterating `S` yields the groups, by increasing primary key and number;
  they are constructed one at a time,
- `collect(S)` returns them as a vector,
- `first(S)` returns the first group, `isempty(S)` tells whether there is one,
- `length(S)` returns the number of groups,
- `keys(S)` iterates over the identifiers of the groups.

`length` and `keys` construct groups only if the library cannot decide a
filter from precomputed data.

# Examples
```jldoctest
julia> L = small_groups_library();

julia> S = find(L, order => 16, !is_abelian)
Selection of small groups: order => 16, !is_abelian

julia> length(S)
9

julia> first(S)
Pc group of order 16

julia> collect(keys(S))
9-element Vector{Tuple{ZZRingElem, ZZRingElem}}:
 (16, 3)
 (16, 4)
 (16, 6)
 (16, 7)
 (16, 8)
 (16, 9)
 (16, 11)
 (16, 12)
 (16, 13)

julia> [describe(G) for G in find(L, 8)]
5-element Vector{String}:
 "C8"
 "C4 x C2"
 "D8"
 "Q8"
 "C2 x C2 x C2"

julia> length(find(L, order => 512, is_abelian))
30
```
"""
function find(L::GroupLibrary, filters...)
  key = _primary_key(L)
  K = _key_type(L)
  keyvals = nothing
  rest = []
  for f in _expand_key_shorthand(filters, key)
    if f isa Pair && _is_primary_key(L, f[1])
      vals = f[2] isa IntegerUnion ? [f[2]] : f[2]
      @req vals isa AbstractVector{<:IntegerUnion} "bad argument $(f[2]) for function $key"
      keyvals = keyvals === nothing ? collect(vals) : intersect(keyvals, vals)
    else
      push!(rest, f)
    end
  end
  @req keyvals !== nothing "must specify a filter for $key"

  keyvals = sort!(unique!(K[k for k in keyvals if k >= 1]))
  rest = Tuple(rest)
  return GroupLibrarySelection(L, keyvals, rest, _translate(L, rest))
end

@doc raw"""
    identify(L::GroupLibrary, G)

Return the identifier `(k, i)` of the group `G` in the library `L`,
that is, `L[k, i]` and `G` are equivalent in the sense of `L`
(isomorphic, or permutation isomorphic).
An exception is thrown if `G` does not belong to `L` or
if `has_identification(L, k)` is `false`,
where `k` is the value of the primary key of `L` for `G`.

# Examples
```jldoctest
julia> identify(small_groups_library(), symmetric_group(4))
(24, 12)

julia> identify(primitive_groups_library(), symmetric_group(4))
(4, 2)
```
"""
function identify end

@doc raw"""
    has_groups(L::GroupLibrary, k::IntegerUnion)

Return whether the groups in `L` with primary key `k` are available.

# Examples
```jldoctest
julia> has_groups(small_groups_library(), 768)
true

julia> has_groups(small_groups_library(), 1024)
false
```
"""
function has_groups end

@doc raw"""
    has_number_of_groups(L::GroupLibrary, k::IntegerUnion)

Return whether the number of groups in `L` with primary key `k` is available,
that is, whether `length(find(L, k))` works.

# Examples
```jldoctest
julia> has_number_of_groups(small_groups_library(), 1024)
true

julia> length(find(small_groups_library(), 1024))
49487367289
```
"""
function has_number_of_groups end

@doc raw"""
    has_identification(L::GroupLibrary, k::IntegerUnion)

Return whether [`identify`](@ref) is available for the groups in `L`
with primary key `k`.

# Examples
```jldoctest
julia> has_identification(small_groups_library(), 256)
true

julia> has_identification(small_groups_library(), 512)
false
```
"""
function has_identification end

function Base.getindex(L::GroupLibrary, id::Tuple{IntegerUnion,IntegerUnion})
  return L[id[1], id[2]]
end

Base.show(io::IO, L::GroupLibrary) = print(io, "Library of ", _name(L))

################################################################################
#
#  Selections
#
################################################################################

Base.IteratorSize(::Type{<:GroupLibrarySelection}) = Base.SizeUnknown()

function Base.eltype(::Type{<:GroupLibrarySelection{<:GroupLibrary{T}}}) where {T}
  return T
end

function Base.iterate(S::GroupLibrarySelection)
  groups = Iterators.flatten((G for (_, G) in _pairs(S, k)) for k in S.keys)
  return _advance(groups, iterate(groups))
end

function Base.iterate(S::GroupLibrarySelection, (groups, state))
  return _advance(groups, iterate(groups, state))
end

_advance(groups, next) = next === nothing ? nothing : (next[1], (groups, next[2]))

function Base.length(S::GroupLibrarySelection)
  return sum(k -> Int(_count(S, k)), S.keys; init=0)
end

Base.keys(S::GroupLibrarySelection) = Iterators.flatten(_ids(S, k) for k in S.keys)

Base.isempty(S::GroupLibrarySelection) = all(k -> _first(S, k) === nothing, S.keys)

function Base.first(S::GroupLibrarySelection)
  for k in S.keys
    G = _first(S, k)
    G === nothing || return G
  end
  throw(ArgumentError("collection must be non-empty"))
end

_show_filter(io::IO, f) = print(io, f)
_show_filter(io::IO, f::ComposedFunction{typeof(!)}) = print(io, "!", f.inner)

function _show_filter(io::IO, f::Pair)
  _show_filter(io, f[1])
  print(io, " => ", f[2])
end

function Base.show(io::IO, S::GroupLibrarySelection)
  lib = S.library
  print(io, "Selection of ", _name(lib), ": ", _primary_key(lib), " => ")
  print(io, length(S.keys) == 1 ? S.keys[1] : S.keys)
  for f in S.rest
    print(io, ", ")
    _show_filter(io, f)
  end
end
