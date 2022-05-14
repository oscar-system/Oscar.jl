export
    all_small_groups,
    has_number_small_groups,
    has_small_group_identification,
    has_small_groups,
    number_small_groups,
    small_group,
    small_group_identification


###################################################################
# Small groups
###################################################################

"""
    has_number_small_groups(n::IntegerUnion)

Return `true` if the number of groups of order `n` is known, otherwise `false`.

# Examples
```jldoctest
julia> has_number_small_groups(1024)
true

julia> has_number_small_groups(2048)
false
```
"""
function has_number_small_groups(n::IntegerUnion)
    n >= 1 || throw(ArgumentError("group order must be positive, not $n"))
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
    n >= 1 || throw(ArgumentError("group order must be positive, not $n"))
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
    n >= 1 || throw(ArgumentError("group order must be positive, not $n"))
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
<pc group of size 60 with 4 generators>

julia> small_group(60, 5)
Group([ (1,2,3,4,5), (1,2,3) ])

julia> small_group(PcGroup, 60, 4)
<pc group of size 60 with 4 generators>
```
"""
function small_group(::Type{T}, n::IntegerUnion, m::IntegerUnion) where T
  G = _small_group(n, m)
  return T(G)
end

function small_group(n::IntegerUnion, m::IntegerUnion)
  G = _small_group(n, m)
  T = _get_type(G)
  return T(G)
end

function _small_group(n::IntegerUnion, m::IntegerUnion)
  N = number_small_groups(n)
  m <= N || throw(ArgumentError("There are only $N groups of order $n, up to isomorphism."))
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
ERROR: identification is not available for groups of order 2432902008176640000
```
"""
function small_group_identification(G::GAPGroup)
   isfinite(G) || error("group is not finite")
   has_small_group_identification(order(G)) || error("identification is not available for groups of order $(order(G))")
   res = GAP.Globals.IdGroup(G.X)
   return Tuple{fmpz,fmpz}(res)
end


"""
    number_small_groups(n::IntegerUnion)

Return the number of groups of order `n`, up to isomorphism.

# Examples
```jldoctest
julia> number_small_groups(8)
5

julia> number_small_groups(4096)
ERROR: the number of groups of order 4096 is not available

julia> number_small_groups(next_prime(fmpz(2)^64))
1
```
"""
function number_small_groups(n::IntegerUnion)
  has_number_small_groups(n) || error("the number of groups of order $n is not available")
  return fmpz(GAP.Globals.NumberSmallGroups(GAP.Obj(n))::GapInt)
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
- `is_abelian`
- `isalmostsimple`
- `is_cyclic`
- `is_nilpotent`
- `isperfect`
- `is_quasisimple`
- `is_simple`
- `is_sporadic_simple`
- `issolvable`
- `issupersolvable`
- `number_conjugacy_classes`
- `order`

The type of the returned groups is `PcGroup` if the group is solvable, `PermGroup` otherwise.

# Examples

List all abelian non-cyclic groups of order 12:
```jldoctest
julia> all_small_groups(12, !is_cyclic, is_abelian)
1-element Vector{PcGroup}:
 <pc group of size 12 with 3 generators>
```

List groups of order 1 to 10 which are not abelian:
```jldoctest
julia> all_small_groups(1:10, !is_abelian)
4-element Vector{PcGroup}:
 <pc group of size 6 with 2 generators>
 <pc group of size 8 with 3 generators>
 <pc group of size 8 with 3 generators>
 <pc group of size 10 with 2 generators>
```
"""
function all_small_groups(L...)
   !isempty(L) || throw(ArgumentError("must specify at least one filter"))
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (order => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L)
   K = GAP.Globals.AllSmallGroups(gapargs...)

   # TODO: perhaps add an option so that ids are returned instead of groups,
   # by calling GAP.Globals.IdsOfAllSmallGroups
   return [_get_type(x)(x) for x in K]
end

#T problem:

#T all_small_groups( 60, is_simple ) -> Array{PermGroup,1}
#T all_small_groups( 60 )  -> Array{Oscar.GAPGroup,1}
#T all_small_groups( 59 )  -> Array{PcGroup,1}

#T Do we want this???
