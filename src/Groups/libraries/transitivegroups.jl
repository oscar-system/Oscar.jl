export
    all_transitive_groups,
    number_transitive_groups,
    transitive_groups_available,
    transitive_group,
    transitive_identification


"""
    transitive_groups_available(deg::Int)

Return whether the transitive groups groups of degree `deg` are available for
use. This function should be used to test for the scope of the library
available.
"""
transitive_groups_available(deg::Int) = GAP.Globals.TransitiveGroupsAvailable(deg)::Bool


"""
    number_transitive_groups(deg::Int)

Return the number of transitive groups of degree `deg` stored in the library
of transitive groups. The function returns `missing` if `deg` is beyond the
range of the library.
"""
function number_transitive_groups(deg::Int)
   res = GAP.Globals.NrTransitiveGroups(deg)
   res isa Int && return res
   return missing
end

"""
    transitive_group(deg::Int, i::Int)

Return the `i`-th group in the catalogue of transitive groups over the set
`{1, ..., deg}` in GAP's Transitive Groups Library.
The output is a group of type `PermGroup`.
"""
function transitive_group(deg::Int, n::Int)
   transitive_groups_available(deg) || error("Transitive groups of degree $(deg) are not available")
   N = number_transitive_groups(deg)
   @assert n <= N "There are only $N transitive groups of degree $deg, up to permutation isomorphism."

   return PermGroup(GAP.Globals.TransitiveGroup(deg,n), deg)
end

"""
    transitive_identification(G::PermGroup)

Return `m` such that `G` is permutation isomorphic with
`transitive_group(d, m)`, where `G` is transitive on `d` points.

If `G` is not transitive, or if the transitive groups of degree `d`
are not available, an exception is raised.

# Examples
```jldoctest
julia> G = symmetric_group(7);  m = transitive_identification(G)
7

julia> order(transitive_group(7, m)) == order(G)
true

julia> S = sub(G, [gap_perm([1, 3, 4, 5, 2])])[1]
Group([ (2,3,4,5) ])

julia> istransitive(S)
false

julia> istransitive(S, moved_points(S))
true

julia> m = transitive_identification(S)
1

julia> order(transitive_group(4, m)) == order(S)
true

julia> transitive_identification(symmetric_group(32))
ERROR: Transitive groups of degree 32 are not available

julia> S = sub(G, [gap_perm([1,3,4,5,2,7,6])])[1];

julia> istransitive(S, moved_points(S))
false

julia> transitive_identification(S)
ERROR: group is not transitive
```
"""
function transitive_identification(G::PermGroup)
  deg = degree(G)
  moved = moved_points(G)
  istransitive(G, moved) || error("group is not transitive")
  transitive_groups_available(deg) || error("Transitive groups of degree $(deg) are not available")
  return GAP.Globals.TransitiveIdentification(G.X)::Int
end

"""
    all_transitive_groups(L...)

Return the list of all transitive groups (up to permutation isomorphism)
satisfying the conditions describe by the arguments. These conditions
may be of one of the following forms:

- `func => intval` selects groups for which the function `func` returns `intval`
- `func => list` selects groups for which the function `func` returns any element inside `list`
- `predicate` selects groups for which the function `predicate` returns `true`
- `!predicate` selects groups for which the function `predicate` returns `false`

The type of the returned groups is `PermGroup`.

TODO: specify which predicates are supported

# Examples
```jldoctest
julia> all_transitive_groups(degree => 4, isabelian)
2-element Vector{PermGroup}:
 C(4) = 4
 E(4) = 2[x]2
```
returns the list of all abelian transitive groups acting on a set of order 4.

The type of the groups is `PermGroup`.
"""
function all_transitive_groups(L...)
   gapargs = translate_group_library_args(L; permgroups=true)
   K = GAP.Globals.AllTransitiveGroups(gapargs...)
   return [PermGroup(x) for x in K]
end

