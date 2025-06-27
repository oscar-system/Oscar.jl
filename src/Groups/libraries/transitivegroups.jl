"""
    has_number_of_transitive_groups(deg::Int)

Return whether the number of transitive groups of degree `deg` is available for
use via `number_of_transitive_groups`.

# Examples
```jldoctest
julia> has_number_of_transitive_groups(30)
true

julia> has_number_of_transitive_groups(64)
false
```
"""
has_number_of_transitive_groups(deg::Int) = has_transitive_groups(deg)

"""
    has_transitive_group_identification(deg::Int)

Return whether identification of transitive groups of degree `deg` is available
via `transitive_group_identification`.

# Examples
```jldoctest
julia> has_transitive_group_identification(30)
true

julia> has_transitive_group_identification(64)
false
```
"""
has_transitive_group_identification(deg::Int) = has_transitive_groups(deg)

"""
    has_transitive_groups(deg::Int)

Return whether the transitive groups of degree `deg` are available for
use. This function should be used to test for the scope of the library
available.

# Examples
```jldoctest
julia> has_transitive_groups(30)
true

julia> has_transitive_groups(64)
false
```
"""
function has_transitive_groups(deg::Int)
  @req deg >= 1 "degree must be positive, not $deg"
  return deg == 1 || GAP.Globals.TransitiveGroupsAvailable(deg)::Bool
end


"""
    number_of_transitive_groups(deg::Int)

Return the number of transitive groups of degree `deg`,
up to permutation isomorphism.

# Examples
```jldoctest
julia> number_of_transitive_groups(30)
5712

julia> number_of_transitive_groups(64)
ERROR: ArgumentError: the number of transitive groups of degree 64 is not available
[...]
```
"""
function number_of_transitive_groups(deg::Int)
  @req has_number_of_transitive_groups(deg) "the number of transitive groups of degree $(deg) is not available"
  return GAP.Globals.NrTransitiveGroups(deg)::Int
end

"""
    transitive_group(deg::Int, i::Int)

Return the `i`-th group in the catalogue of transitive groups over the set
`{1, ..., deg}` in GAP's Transitive Groups Library.
The output is a group of type `PermGroup`.

# Examples
```jldoctest
julia> transitive_group(5,4)
Alternating group of degree 5 with 2 generators
  (1,2,3,4,5)
  (3,4,5)

julia> transitive_group(5,6)
ERROR: ArgumentError: there are only 5 transitive groups of degree 5, not 6
[...]
```
"""
function transitive_group(deg::Int, i::Int)
  @req has_transitive_groups(deg) "transitive groups of degree $deg are not available"
  N = number_of_transitive_groups(deg)
  @req i <= N "there are only $N transitive groups of degree $deg, not $i"
  return PermGroup(GAP.Globals.TransitiveGroup(deg,i), deg)
end

"""
    transitive_group_identification(G::PermGroup)

Return a pair `(d,n)` such that `G` is permutation isomorphic with
`transitive_group(d,n)`, where `G` acts transitively on `d` points.

If `G` is not transitive on its moved points, or if the transitive groups of
degree `d` are not available, an exception is thrown.

# Examples
```jldoctest
julia> G = symmetric_group(7);  m = transitive_group_identification(G)
(7, 7)

julia> order(transitive_group(m...)) == order(G)
true

julia> S = sub(G, [perm([1, 3, 4, 5, 2])])[1]
Permutation group of degree 7 with 1 generator
  (2,3,4,5)

julia> is_transitive(S)
false

julia> is_transitive(S, moved_points(S))
true

julia> m = transitive_group_identification(S)
(4, 1)

julia> order(transitive_group(m...)) == order(S)
true

julia> transitive_group_identification(symmetric_group(64))
ERROR: ArgumentError: identification of transitive groups of degree 64 are not available
[...]

julia> S = sub(G, [perm([1,3,4,5,2,7,6])])[1];

julia> transitive_group_identification(S)
ERROR: ArgumentError: group is not transitive on its moved points
[...]
```
"""
function transitive_group_identification(G::PermGroup)
  moved = moved_points(G)
  @req is_transitive(G, moved) "group is not transitive on its moved points"
  deg = length(moved)
  @req has_transitive_groups(deg) "identification of transitive groups of degree $(deg) are not available"
  res = GAP.Globals.TransitiveIdentification(GapObj(G))::Int
  return deg, res
end

"""
    all_transitive_groups(L...)

Return the list of all transitive groups (up to permutation isomorphism)
satisfying the conditions described by the arguments. These conditions
may be of one of the following forms:

- `func => intval` selects groups for which the function `func` returns `intval`
- `func => list` selects groups for which the function `func` returns any element inside `list`
- `func` selects groups for which the function `func` returns `true`
- `!func` selects groups for which the function `func` returns `false`

As a special case, the first argument may also be one of the following:
- `intval` selects groups whose degree equals `intval`; this is equivalent to `degree => intval`
- `intlist` selects groups whose degree is in `intlist`; this is equivalent to `degree => intlist`

The following functions are currently supported as values for `func`:
- `degree`
- `is_abelian`
- `is_almost_simple`
- `is_cyclic`
- `is_nilpotent`
- `is_perfect`
- `is_primitive`
- `is_quasisimple`
- `is_simple`
- `is_sporadic_simple`
- `is_solvable`
- `is_supersolvable`
- `is_transitive`
- `number_of_conjugacy_classes`
- `number_of_moved_points`
- `order`
- `transitivity`

The type of the returned groups is `PermGroup`.

# Examples
```jldoctest
julia> all_transitive_groups(4)
5-element Vector{PermGroup}:
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 4
 Alternating group of degree 4
 Symmetric group of degree 4

julia> all_transitive_groups(degree => 3:5, is_abelian)
4-element Vector{PermGroup}:
 Alternating group of degree 3
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 5
```
"""
function all_transitive_groups(L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (degree => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L; filter_attrs = _permgroup_filter_attrs)
   K = GAP.Globals.AllTransitiveGroups(gapargs...)
   return [PermGroup(x) for x in K]
end

# TODO: turn this into an iterator, possibly using PrimitiveGroupsIterator
