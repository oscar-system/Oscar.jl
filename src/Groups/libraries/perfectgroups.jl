export
    has_number_perfect_groups,
    has_perfect_group_identification,
    has_perfect_groups,
    number_perfect_groups,
    orders_perfect_groups,
    perfect_group_identification,
    perfect_group


###################################################################
# Perfect groups
###################################################################

@doc Markdown.doc"""
    orders_perfect_groups()

Returns a sorted vector of all numbers to $2 \cdot 10^6$ that occur as orders
of perfect groups.

# Examples
```jldoctest
julia> orders_perfect_groups()[1:10]
10-element Vector{Int64}:
   1
  60
 120
 168
 336
 360
 504
 660
 720
 960
```
"""
function orders_perfect_groups()
    return Vector{Int}(GAP.Globals.SizesPerfectGroups())
end

@doc Markdown.doc"""
    has_number_perfect_groups(n::Int)

Return `true` if the number of perfect groups of order `n` are available
via `number_perfect_groups`, otherwise `false`.

Currently the number of perfect groups is available up to order $2 \cdot 10^6$.

# Examples
```jldoctest
julia> has_number_perfect_groups(7200)
true

julia> has_number_perfect_groups(2*10^6+1)
false
```
"""
has_number_perfect_groups(n::Int) = has_perfect_groups(n) # for now just an alias

@doc Markdown.doc"""
    has_perfect_group_identification(n::Int)

Return `true` if identification is supported for the perfect
groups of order `n` via `perfect_group_identification`, otherwise `false`.

Currently identification is available for all perfect groups up to order $2 \cdot 10^6$.

# Examples
```jldoctest
julia> has_perfect_group_identification(7200)
true

julia> has_perfect_group_identification(2*10^6+1)
false
```
"""
has_perfect_group_identification(n::Int) = has_perfect_groups(n) # for now just an alias

@doc Markdown.doc"""
    has_perfect_groups(deg::Int)

Return `true` if the perfect groups of order `n` are available
via `perfect_group` and `all_perfect_groups`, otherwise `false`.

Currently all perfect groups up to order $2 \cdot 10^6$ are available.

# Examples
```jldoctest
julia> has_perfect_groups(7200)
true

julia> has_perfect_groups(2*10^6+1)
false
```
"""
function has_perfect_groups(n::Int)
  n >= 1 || throw(ArgumentError("order must be positive, not $n"))
  return n <= 2_000_000
end

"""
    perfect_group(::Type{T} = PermGroup, n::IntegerUnion, k::IntegerUnion)

Return the `k`-th group of order `n` and type `T` in the catalogue of
perfect groups in GAP's Perfect Groups Library.
The type `T` can be either `PermGroup` or `FPGroup`.

# Examples
```jldoctest
julia> perfect_group(60, 1)
A5

julia> gens(ans)
2-element Vector{PermGroupElem}:
 (1,2)(4,5)
 (2,3,4)

julia> perfect_group(FPGroup, 60, 1)
A5

julia> gens(ans)
2-element Vector{FPGroupElem}:
 a
 b

julia> perfect_group(60, 2)
ERROR: ArgumentError: there are only 1 perfect groups of order 60
```
"""
function perfect_group(::Type{T}, n::IntegerUnion, k::IntegerUnion) where T <: GAPGroup
   n >= 1 || throw(ArgumentError("group order must be positive, not $n"))
   k >= 1 || throw(ArgumentError("group index must be positive, not $k"))
   N = number_perfect_groups(n)
   k <= N || throw(ArgumentError("there are only $N perfect groups of order $n"))
   if T == PermGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsPermGroup, GAP.Obj(n), GAP.Obj(k)))
   elseif T == FPGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsSubgroupFpGroup, GAP.Obj(n), GAP.Obj(k)))
   else
      throw(ArgumentError("unsupported type $T"))
   end

   return G
end

perfect_group(n::IntegerUnion, m::IntegerUnion) = perfect_group(PermGroup, n, m)

"""
    perfect_group_identification(G::GAPGroup)

Return `(n, m)` such that `G` is isomorphic with `perfect_group(n, m)`.
If `G` is not perfect, an exception is thrown.

# Examples
```jldoctest
julia> perfect_group_identification(alternating_group(5))
(60, 1)

julia> perfect_group_identification(SL(2,7))
(336, 1)
```
"""
function perfect_group_identification(G::GAPGroup)
   is_perfect(G) || error("group is not perfect")
   res = GAP.Globals.PerfectIdentification(G.X)
   res !== GAP.Globals.fail || error("identification is not available for groups of order $(order(G))")
   return Tuple{Int,Int}(res)
end

"""
    number_perfect_groups(n::IntegerUnion)

Return the number of perfect groups of order `n`, up to isomorphism.

# Examples
```jldoctest
julia> number_perfect_groups(60)
1

julia> number_perfect_groups(1966080)
7344
```
"""
function number_perfect_groups(n::IntegerUnion)
   n >= 1 || throw(ArgumentError("group order must be positive, not $n"))
   res = GAP.Globals.NumberPerfectGroups(GAP.Obj(n))
   res !== GAP.Globals.fail || error("the number of perfect groups of order $n is not available")
   return res::Int
end

# TODO: add all_perfect_groups() iterator
