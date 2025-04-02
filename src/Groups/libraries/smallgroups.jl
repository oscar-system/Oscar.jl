###################################################################
# Small groups
###################################################################

"""
    has_number_of_small_groups(n::IntegerUnion)

Return `true` if the number of groups of order `n` is known, otherwise `false`.

# Examples
```jldoctest
julia> has_number_of_small_groups(1024)
true

julia> has_number_of_small_groups(2048)
false
```
"""
function has_number_of_small_groups(n::IntegerUnion)
    @req n >= 1 "group order must be positive, not $n"
    return GAP.Globals.NumberSmallGroupsAvailable(GAP.Obj(n))
end

"""
    has_small_groups(n::IntegerUnion)

Return `true` if the groups of order `n` are available via `small_group`
and `all_small_groups`, otherwise `false`.

# Examples
```jldoctest
julia> has_small_groups(512)
true

julia> has_small_groups(1024)
false
```
"""
function has_small_groups(n::IntegerUnion)
    @req n >= 1 "group order must be positive, not $n"
    return GAP.Globals.SmallGroupsAvailable(GAP.Obj(n))
end

"""
    has_small_group_identification(n::IntegerUnion)

Return `true` if identification for groups of order `n` is available via
`small_group_identification`, otherwise `false`.

# Examples
```jldoctest
julia> has_small_group_identification(256)
true

julia> has_small_group_identification(512)
false
```
"""
function has_small_group_identification(n::IntegerUnion)
    @req n >= 1 "group order must be positive, not $n"
    return GAP.Globals.IdGroupsAvailable(GAP.Obj(n))
end

"""
    small_group(::Type{T}, n::IntegerUnion, i::IntegerUnion) where T
    small_group(n::IntegerUnion, i::IntegerUnion)

Return the `i`-th group of order `n` in the Small Groups Library. If a type
`T` is specified then an attempt is made to return the result with that type.
If `T` is omitted then the resulting group will have type `PcGroup` if it is
solvable, otherwise it will be of type `PermGroup`.

# Examples
```jldoctest
julia> small_group(60, 4)
Pc group of order 60

julia> small_group(60, 5)
Permutation group of degree 5 and order 60

julia> small_group(PcGroup, 60, 4)
Pc group of order 60
```
"""
function small_group(::Type{T}, n::IntegerUnion, m::IntegerUnion) where T
  G = small_group(n, m)
  return T(G)
end

function small_group(n::IntegerUnion, m::IntegerUnion)
  G = _small_group(n, m)
  return _oscar_group(G)
end

function _small_group(n::IntegerUnion, m::IntegerUnion)
  N = number_of_small_groups(n)
  @req m <= N "There are only $N groups of order $n, up to isomorphism."
  return GAP.Globals.SmallGroup(GAP.Obj(n), GAP.Obj(m))
end


"""
    small_group_identification(G::Group)

Return a pair of integer `(n, m)`, where `G` is isomorphic with `small_group(n, m)`.

# Examples
```jldoctest
julia> small_group_identification(alternating_group(5))
(60, 5)

julia> small_group_identification(symmetric_group(20))
ERROR: ArgumentError: identification is not available for groups of order 2432902008176640000
[...]
```
"""
function small_group_identification(G::GAPGroup)
   @req is_finite(G) "group is not finite"
   @req has_small_group_identification(order(G)) "identification is not available for groups of order $(order(G))"
   res = GAP.Globals.IdGroup(GapObj(G))
   return Tuple{ZZRingElem,ZZRingElem}(res)
end


"""
    number_of_small_groups(n::IntegerUnion)

Return the number of groups of order `n`, up to isomorphism.

# Examples
```jldoctest
julia> number_of_small_groups(8)
5

julia> number_of_small_groups(4096)
ERROR: ArgumentError: the number of groups of order 4096 is not available
[...]

julia> number_of_small_groups(next_prime(ZZRingElem(2)^64))
1
```
"""
function number_of_small_groups(n::IntegerUnion)
  @req has_number_of_small_groups(n) "the number of groups of order $n is not available"
  return ZZRingElem(GAP.Globals.NumberSmallGroups(GAP.Obj(n))::GapInt)
end


"""
    all_small_groups(L...)

Return the list of all groups (up to isomorphism)
satisfying the conditions described by the arguments. These conditions
may be of one of the following forms:

- `func => intval` selects groups for which the function `func` returns `intval`
- `func => list` selects groups for which the function `func` returns any element inside `list`
- `func` selects groups for which the function `func` returns `true`
- `!func` selects groups for which the function `func` returns `false`

As a special case, the first argument may also be one of the following:
- `intval` selects groups whose order equals `intval`; this is equivalent to `order => intval`
- `intlist` selects groups whose order is in `intlist`; this is equivalent to `order => intlist`

Note that at least one of the conditions must impose a limit on the group
order, otherwise an exception is thrown.

The following functions are currently supported as values for `func`:
- `exponent`
- `is_abelian`
- `is_almost_simple`
- `is_cyclic`
- `is_nilpotent`
- `is_perfect`
- `is_quasisimple`
- `is_simple`
- `is_sporadic_simple`
- `is_solvable`
- `is_supersolvable`
- `number_of_conjugacy_classes`
- `order`

The type of the returned groups is `PcGroup` if the group is solvable, `PermGroup` otherwise.

# Examples

List all abelian non-cyclic groups of order 12:
```jldoctest
julia> all_small_groups(12, !is_cyclic, is_abelian)
1-element Vector{PcGroup}:
 Pc group of order 12
```

List groups of order 1 to 10 which are not abelian:
```jldoctest
julia> all_small_groups(1:10, !is_abelian)
4-element Vector{PcGroup}:
 Pc group of order 6
 Pc group of order 8
 Pc group of order 8
 Pc group of order 10
```
"""
function all_small_groups(L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (order => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L)
   K = GAP.Globals.AllSmallGroups(gapargs...)

   # TODO: perhaps add an option so that ids are returned instead of groups,
   # by calling GAP.Globals.IdsOfAllSmallGroups
   return [_oscar_group(x) for x in K]
end

#T problem:

#T all_small_groups( 60, is_simple ) -> Array{PermGroup,1}
#T all_small_groups( 60 )  -> Array{Oscar.GAPGroup,1}
#T all_small_groups( 59 )  -> Array{PcGroup,1}

#T Do we want this???
