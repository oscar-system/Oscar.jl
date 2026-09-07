@doc raw"""
The functions in this module are wrappers for the GAP library of
finite perfect groups which provides, up to isomorphism,
a list of all perfect groups whose sizes are less than $2 \cdot 10^6$.
The groups of most orders up to $10^6$ have been
enumerated by Derek Holt and Wilhelm Plesken, see [HP89](@cite). For orders
86016, 368640, or 737280 this work only counted the groups (but did not
explicitly list them), the groups of orders 61440, 122880, 172032,
245760, 344064, 491520, 688128, or 983040 were omitted.

Several additional groups omitted from the book [HP89](@cite) have also
been included. Two groups -- one of order 450000 with a factor group of
type $A_6$ and the one of order 962280 -- were found in 2005 by Jack Schmidt.
Two groups of order 243000 and one each of orders 729000, 871200, 878460
were found in 2020 by Alexander Hulpke.

The perfect groups of size less than $2 \cdot 10^6$ which had not been
classified in the work of Holt and Plesken have been enumerated by Alexander
Hulpke, see [Hul22](@cite).
They are stored directly and provide less construction information
in their names.

As all groups are stored by presentations, a permutation representation
is obtained by coset enumeration. Note that some of the library groups do
not have a faithful permutation representation of small degree.
Computations in these groups may be rather time consuming.
"""
module PerfectGroups

using Oscar
import Oscar.GAPGroup
import Oscar.IntegerUnion

@doc raw"""
    PerfectGroups.has_count(n::Int)

Return whether the number of perfect groups of order `n` is available
via `PerfectGroups.count(n)`.

Currently the number of perfect groups is available up to order $2 \cdot 10^6$.

# Examples
```jldoctest
julia> PerfectGroups.has_count(7200)
true

julia> PerfectGroups.has_count(2*10^6+1)
false
```
"""
has_count(n::Int) = has(n) # for now just an alias

"""
    PerfectGroups.has_identification(n::Int)

Return whether identification of perfect groups of order `n` is available
via `PerfectGroups.identification(G)`, where `G` is a group of order `n`.

# Examples
```jldoctest
julia> PerfectGroups.has_identification(7200)
true

julia> PerfectGroups.has_identification(2*10^6+1)
false
```
"""
has_identification(n::Int) = has(n) # for now just an alias

@doc raw"""
    PerfectGroups.has(n::Int)

Return whether the perfect groups of order `n` are available for
use.  This function should be used to test for the scope of the library
available.

Currently all perfect groups up to order $2 \cdot 10^6$ are available.

# Examples
```jldoctest
julia> PerfectGroups.has(7200)
true

julia> PerfectGroups.has(2*10^6+1)
false
```
"""
function has(n::Int)
  @req n >= 1 "order must be positive, not $n"
  return n <= 2_000_000
end

"""
    PerfectGroups.count(n::IntegerUnion)

Return the number of perfect groups of order `n`, up to isomorphism.

# Examples
```jldoctest
julia> PerfectGroups.count(60)
1

julia> PerfectGroups.count(1966080)
7344
```
"""
function count(n::IntegerUnion)
   @req n >= 1 "group order must be positive, not $n"
   res = GAP.Globals.NumberPerfectGroups(GAP.Obj(n))
   @req (res !== GAP.Globals.fail) "the number of perfect groups of order $n is not available"
   return res::Int
end

"""
    PerfectGroups.get(::Type{T} = PermGroup, n::IntegerUnion, k::IntegerUnion)

Return the `k`-th group of order `n` and type `T` in the catalogue of
perfect groups in GAP's Perfect Groups Library.
The type `T` can be either `PermGroup` or `FPGroup`.

# Examples
```jldoctest
julia> PerfectGroups.get(60, 1)
Permutation group of degree 5 and order 60

julia> gens(ans)
2-element Vector{PermGroupElem}:
 (1,2)(4,5)
 (2,3,4)

julia> PerfectGroups.get(FPGroup, 60, 1)
Finitely presented group of order 60

julia> gens(ans)
2-element Vector{FPGroupElem}:
 a
 b

julia> PerfectGroups.get(60, 2)
ERROR: ArgumentError: there are only 1 perfect groups of order 60
[...]
```
"""
function get(::Type{T}, n::IntegerUnion, k::IntegerUnion) where T <: GAPGroup
   @req n >= 1 "group order must be positive, not $n"
   @req k >= 1 "group index must be positive, not $k"
   N = count(n)
   @req k <= N "there are only $N perfect groups of order $n"
   if T == PermGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsPermGroup, GAP.Obj(n), GAP.Obj(k)))
   elseif T == FPGroup
      G = T(GAP.Globals.PerfectGroup(GAP.Globals.IsSubgroupFpGroup, GAP.Obj(n), GAP.Obj(k)))
   else
      throw(ArgumentError("unsupported type $T"))
   end

   return G
end

get(n::IntegerUnion, m::IntegerUnion) = get(PermGroup, n, m)

"""
    PerfectGroups.identification(G::GAPGroup)

Return a pair `(n, m)` such that `G` is isomorphic with
`PerfectGroups.get(n, m)`.
If `G` is not perfect, an exception is thrown.

# Examples
```jldoctest
julia> PerfectGroups.identification(alternating_group(5))
(60, 1)

julia> PerfectGroups.identification(SL(2,7))
(336, 1)
```
"""
function identification(G::GAPGroup)
   @req is_perfect(G) "group is not perfect"
   res = GAP.Globals.PerfectIdentification(GapObj(G))
   @req (res !== GAP.Globals.fail) "identification is not available for groups of order $(order(G))"
   return Tuple{Int,Int}(res)
end

"""
    PerfectGroups.get_all(L...)
    PerfectGroups.get_one(L...)

`PerfectGroups.get_all` returns the list of all perfect groups
(up to isomorphism) satisfying the conditions described by the arguments.
`PerfectGroups.get_one` returns one perfect group satisfying
the conditions if such a group exists, and `nothing` otherwise.
These conditions may be of one of the following forms:

- `func => intval` selects groups for which the function `func` returns `intval`
- `func => list` selects groups for which the function `func` returns any element inside `list`
- `func` selects groups for which the function `func` returns `true`
- `!func` selects groups for which the function `func` returns `false`

As a special case, the first argument may also be one of the following:
- `intval` selects groups whose order equals `intval`; this is equivalent to `order => intval`
- `intlist` selects groups whose order is in `intlist`; this is equivalent to `order => intlist`

The following functions are currently supported as values for `func`:
- `is_quasisimple`
- `is_simple`
- `is_sporadic_simple`
- `number_of_conjugacy_classes`
- `order`

The type of the returned group(s) is `PermGroup`.

# Examples
```jldoctest
julia> PerfectGroups.get_all(7200)
2-element Vector{PermGroup}:
 Permutation group of degree 29 and order 7200
 Permutation group of degree 288 and order 7200

julia> PerfectGroups.get_all(order => 1:200, !is_simple)
2-element Vector{PermGroup}:
 Permutation group of degree 1 and order 1
 Permutation group of degree 24 and order 120

julia> PerfectGroups.get_one(order => 1:200, !is_simple)
Permutation group of degree 1 and order 1
```
"""
get_all(L...) = _get_all_or_one(L, true)

get_one(L...) = _get_all_or_one(L, false)
@doc (@doc get_all) get_one

function _get_all_or_one(L, all::Bool)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (order => L[1], L[2:end]...)
   end
   # first get all order restrictions
   ordsL = [x for x in L if x isa Pair && x[1] == order]
   @req !isempty(ordsL) "must restrict the order"
   conds = [x for x in L if !(x isa Pair && x[1] == order)]
   orders = intersect([x[2] for x in ordsL]...)
   @req has_perfect_groups(maximum(orders)) "only orders up to 2 million are supported"
   res = PermGroup[]
   for n in orders, i in 1:number_of_perfect_groups(n)
      G = perfect_group(n, i)
      ok = true
      for c in conds
         if c isa Pair
            val = c[1](G)
            ok = (val == c[2] || val in c[2])
         elseif c isa Function
            ok = c(G)
         else
            throw(ArgumentError("expected a function or a pair, got $arg"))
         end
         ok || break
      end
      !all && ok && return G
      ok && push!(res, G)
   end
   return all ? res : nothing
end

@doc raw"""
    PerfectGroups.orders()

Return a sorted vector of all numbers to $2 \cdot 10^6$ that occur as orders
of perfect groups.

# Examples
```jldoctest
julia> PerfectGroups.orders()[1:10]
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
orders() = Vector{Int}(GAP.Globals.SizesPerfectGroups())

end # module PerfectGroups

export PerfectGroups
