export
    orders_perfect_groups,
    number_perfect_groups,
    perfect_group,
    perfect_group_identification


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
```
"""
function perfect_group(::Type{T}, n::IntegerUnion, k::IntegerUnion) where T <: GAPGroup
   n >= 1 || throw(ArgumentError("group order must be positive, not $n"))
   k >= 1 || throw(ArgumentError("group index must be positive, not $k"))
   N = number_perfect_groups(n)
   k <= N || throw(ArgumentError("There are only $N perfect groups of order $n, up to isomorphism"))
   if T == PermGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsPermGroup, GAP.Obj(n), GAP.Obj(k)))
   elseif T == FPGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsSubgroupFpGroup, GAP.Obj(n), GAP.Obj(k)))
   else
      throw(ArgumentError("Wrong type"))
   end

   return G
end

perfect_group(n::IntegerUnion, m::IntegerUnion) = perfect_group(PermGroup, n, m)

"""
    perfect_group_identification(G::GAPGroup)

Return `(n, m)` such that `G` is isomorphic with `perfect_group(n, m)`.

# Examples
```jldoctest
julia> perfect_group_identification(alternating_group(5))
(60, 1)

julia> perfect_group_identification(SL(2,7))
(336, 1)
```
"""
function perfect_group_identification(G::GAPGroup)
   isperfect(G) || error("group is not perfect")
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
   res !== GAP.Globals.fail || error("the number of perfect groups of order $n is not in the library")
   return res::Int
end
