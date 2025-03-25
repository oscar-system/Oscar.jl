"""
    has_number_of_primitive_groups(deg::Int)

Return `true` if the number of primitive permutation groups of degree `deg` is available
via `number_of_primitive_groups`, otherwise `false`.

Currently the number of primitive permutation groups is available up to degree 4095.

# Examples
```jldoctest
julia> has_number_of_primitive_groups(50)
true

julia> has_number_of_primitive_groups(5000)
false
```
"""
has_number_of_primitive_groups(deg::Int) = has_primitive_groups(deg) # for now just an alias

"""
    has_primitive_group_identification(deg::Int)

Return `true` if identification is supported for the primitive permutation
groups of degree `deg` via `primitive_group_identification`, otherwise `false`.

Currently identification is available for all primitive permutation groups up to degree 4095.

# Examples
```jldoctest
julia> has_primitive_group_identification(50)
true

julia> has_primitive_group_identification(5000)
false
```
"""
has_primitive_group_identification(deg::Int) = has_primitive_groups(deg) # for now just an alias

"""
    has_primitive_groups(deg::Int)

Return `true` if the primitive permutation groups of degree `deg` are available
via `primitive_group` and `all_primitive_groups`, otherwise `false`.

Currently all primitive permutation groups up to degree 4095 are available.

# Examples
```jldoctest
julia> has_primitive_groups(50)
true

julia> has_primitive_groups(5000)
false
```
"""
function has_primitive_groups(deg::Int)
  @req deg >= 1 "degree must be positive, not $deg"
  return GAP.Globals.PrimitiveGroupsAvailable(GAP.Obj(deg))
end

"""
    number_of_primitive_groups(deg::Int)

Return the number of primitive permutation groups of degree `deg`,
up to permutation isomorphism.

# Examples
```jldoctest
julia> number_of_primitive_groups(10)
9

julia> number_of_primitive_groups(4096)
ERROR: ArgumentError: the number of primitive permutation groups of degree 4096 is not available
[...]
```
"""
function number_of_primitive_groups(deg::Int)
  @req has_number_of_primitive_groups(deg) "the number of primitive permutation groups of degree $deg is not available"
  return GAP.Globals.NrPrimitiveGroups(deg)::Int
end

"""
    primitive_group(deg::Int, i::Int)

Return the `i`-th group in the catalogue of primitive permutation groups over the set
`{1, ..., deg}` in GAP's library of primitive permutation groups.
The output is a group of type `PermGroup`.

# Examples
```jldoctest
julia> primitive_group(10,1)
Permutation group of degree 10 and order 60

julia> primitive_group(10,10)
ERROR: ArgumentError: there are only 9 primitive permutation groups of degree 10, not 10
[...]
```
"""
function primitive_group(deg::Int, i::Int)
  @req has_primitive_groups(deg) "primitive permutation groups of degree $deg are not available"
  N = number_of_primitive_groups(deg)
  @req i <= N "there are only $N primitive permutation groups of degree $deg, not $i"
  return PermGroup(GAP.Globals.PrimitiveGroup(deg,i), deg)
end

"""
    primitive_group_identification(G::PermGroup)

Return a pair `(d,n)` such that `G` is permutation isomorphic with
`primitive_group(d,n)`, where `G` acts primitively on `d` points.

If `G` is not primitive on its moved points, or if the primitive permutation groups of
degree `d` are not available, an exception is thrown.

# Examples
```jldoctest
julia> G = symmetric_group(7);  m = primitive_group_identification(G)
(7, 7)

julia> order(primitive_group(m...)) == order(G)
true

julia> S = stabilizer(G, 1)[1]
Permutation group of degree 7 and order 720

julia> is_primitive(S)
false

julia> is_primitive(S, moved_points(S))
true

julia> m = primitive_group_identification(S)
(6, 4)

julia> order(primitive_group(m...)) == order(S)
true

julia> primitive_group_identification(symmetric_group(4096))
ERROR: ArgumentError: identification of primitive permutation groups of degree 4096 is not available
[...]

julia> S = sub(G, [perm([1,3,4,5,2,7,6])])[1];

julia> primitive_group_identification(S)
ERROR: ArgumentError: group is not primitive on its moved points
[...]
```
"""
function primitive_group_identification(G::PermGroup)
  moved = moved_points(G)
  @req is_primitive(G, moved) "group is not primitive on its moved points"
  deg = length(moved)
  @req has_primitive_groups(deg) "identification of primitive permutation groups of degree $(deg) is not available"
  res = GAP.Globals.PrimitiveIdentification(GapObj(G))::Int
  return deg, res
end

"""
    all_primitive_groups(L...)

Return the list of all primitive permutation groups (up to permutation isomorphism)
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
julia> all_primitive_groups(4)
2-element Vector{PermGroup}:
 Alt(4)
 Sym(4)

julia> all_primitive_groups(degree => 3:5, is_abelian)
2-element Vector{PermGroup}:
 Alt(3)
 Permutation group of degree 5 and order 5
```
"""
function all_primitive_groups(L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (degree => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L; filter_attrs = _permgroup_filter_attrs)
   K = GAP.Globals.AllPrimitiveGroups(gapargs...)
   return [PermGroup(x) for x in K]
end

# TODO: turn this into an iterator, possibly using PrimitiveGroupsIterator
