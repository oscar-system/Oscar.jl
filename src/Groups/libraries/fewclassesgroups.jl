"""
    has_number_of_groups_with_class_number(n::IntegerUnion)

Return whether the number of groups with exactly `n` conjugacy classes is
available for use via [`number_of_groups_with_class_number`](@ref).

# Examples
```jldoctest
julia> has_number_of_groups_with_class_number(14)
true

julia> has_number_of_groups_with_class_number(50)
false
```
"""
has_number_of_groups_with_class_number(n::IntegerUnion) = has_groups_with_class_number(n)

"""
    has_groups_with_class_number(n::IntegerUnion)

Return whether the groups with exactly `n` conjugacy classes are available for
use via [`group_with_class_number`](@ref) and
[`all_groups_with_class_number`](@ref).
This function should be used to test for the scope of the library available.

# Examples
```jldoctest
julia> has_groups_with_class_number(14)
true

julia> has_groups_with_class_number(50)
false
```
"""
function has_groups_with_class_number(n::IntegerUnion)
  @req n >= 1 "class number must be positive, not $n"
  return GAP.Globals.SmallClassNrGroupsAvailable(GapObj(n))
end

"""
    number_of_groups_with_class_number(n::IntegerUnion)

Return the number of groups with exactly `n` conjugacy classes,
up to isomorphism.

# Examples
```jldoctest
julia> number_of_groups_with_class_number(8)
21

julia> number_of_groups_with_class_number(50)
ERROR: ArgumentError: the number of groups with 50 classes is not available
[...]
```
"""
function number_of_groups_with_class_number(n::IntegerUnion)
  @req has_number_of_groups_with_class_number(n) "the number of groups with $n classes is not available"
  return GAP.Globals.NrSmallClassNrGroups(GapObj(n))
end

"""
    has_groups_with_class_number_identification(n::IntegerUnion)

Return `true` if identification for groups with `n` conjugacy classes
is available via `group_with_class_number_identification`, otherwise `false`.

# Examples
```jldoctest
julia> has_groups_with_class_number_identification(14)
true

julia> has_groups_with_class_number_identification(50)
false
```
"""
function has_groups_with_class_number_identification(n::IntegerUnion)
  return has_groups_with_class_number(n)
end

"""
    group_with_class_number(::Type{T}, n::IntegerUnion, i::IntegerUnion) where T
    group_with_class_number(n::IntegerUnion, i::IntegerUnion)

Return the `i`-th group in the list of groups with exactly `n`
conjugacy classes.
If a type `T` is specified then an attempt is made to return the result
with that type.
If `T` is omitted then the resulting group will have type `PcGroup`
if it is a solvable group that belongs to the library of small groups,
otherwise it will be of type `PermGroup`.

# Examples
```jldoctest
julia> describe(group_with_class_number(5, 4))
"D14"

julia> describe(group_with_class_number(5, 8))
"A5"

julia> group_with_class_number(5, 12)
ERROR: ArgumentError: there are only 8 groups with 5 classes, up to isomorphism, not 12
[...]
```
"""
function group_with_class_number(::Type{T}, n::IntegerUnion, i::IntegerUnion) where T
  return T(group_with_class_number(n, i))
end

function group_with_class_number(n::IntegerUnion, i::IntegerUnion)
  @req has_groups_with_class_number(n) "groups with $n classes are not available"
  @req i >= 1 "index must be positive, not $i"
  N = number_of_groups_with_class_number(n)
  @req i <= N "there are only $N groups with $n classes, up to isomorphism, not $i"
  return _oscar_group(GAP.Globals.SmallClassNrGroup(GapObj(n), GapObj(i)))
end

"""
    group_with_class_number_identification(G::GapGroup)

Return a pair of integer `(n, m)`, where `G` is isomorphic with
`group_with_class_number(n, m)`.

# Examples
```jldoctest
julia> group_with_class_number_identification(alternating_group(4))
(4, 4)

julia> group_with_class_number_identification(symmetric_group(20))
ERROR: ArgumentError: identification is not available for groups with 627 conjugacy classes
[...]
```
"""
function group_with_class_number_identification(G::GAPGroup)
   @req is_finite(G) "group is not finite"
   k = number_of_conjugacy_classes(G)
   @req has_groups_with_class_number_identification(k) "identification is not available for groups with $k conjugacy classes"
   res = GAP.Globals.IdClassNr(GapObj(G))
   return Tuple{ZZRingElem,ZZRingElem}(res)
end

"""
    all_groups_with_class_number(L...)

Return the vector of all groups (up to isomorphism) satisfying the conditions
described by the arguments.  These conditions and the format are the same
as for [`all_small_groups`](@ref),
except that the special case where the first argument is an integer or a
list of integers has the following meaning.

- `intval` selects groups with exactly `intval` conjugacy classes;
  this is equivalent to `number_of_conjugacy_classes => intval`

- `intlist` selects groups whose number of conjugacy classes is in `intlist`;
  this is equivalent to `number_of_conjugacy_classes => intlist`

Note that at least one of the conditions must impose a limit on the
number of conjugacy classes, otherwise an exception is thrown.

The type of the returned groups is `PcGroup` if the group is solvable,
`PermGroup` otherwise.

# Examples

List all groups with exactly 4 conjugacy classes:

```jldoctest
julia> map(describe, all_groups_with_class_number(4))
4-element Vector{String}:
 "C4"
 "C2 x C2"
 "D10"
 "A4"
```

List all groups with up to 4 conjugacy classes that are not abelian:

```jldoctest
julia> map(describe, all_groups_with_class_number(1:4, !is_abelian))
3-element Vector{String}:
 "S3"
 "D10"
 "A4"
```
"""
function all_groups_with_class_number(L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (number_of_conjugacy_classes => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L)
   K = GAP.Globals.AllSmallClassNrGroups(gapargs...)

   return [_oscar_group(x) for x in K]
end
