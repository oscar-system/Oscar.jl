"""
The functions in this module are wrappers for
the GAP library of transitive permutation groups up to degree 48,
via the GAP package `TransGrp` [Hul23](@cite).

(The groups of degrees 32 and 48 are currently not automatically
available in OSCAR,
one has to install additional data in order to access them.)

The arrangement and the names of the groups of degree up to 15
is the same as given in [CHM98](@cite).
With the exception of the symmetric and alternating groups
(which are represented as `symmetric_group` and `alternating_group`)
the generators for these groups also conform to this paper
with the only difference that 0 (which is not permitted in GAP
for permutations to act on) is always replaced by the degree.

The arrangement for all degrees is intended to be equal to
the arrangement within the systems GAP and Magma,
thus it should be safe to refer to particular (classes of) groups
by their index numbers.
"""
module TransitiveGroups

using Oscar
import Oscar.IntegerUnion
import Oscar._permgroup_filter_attrs
import Oscar.translate_group_library_args

"""
    TransitiveGroups.has_count(deg::Int)

Return whether the number of transitive groups of degree `deg` is available for
use via `TransitiveGroups.count(deg)`.

# Examples
```jldoctest
julia> TransitiveGroups.has_count(30)
true

julia> TransitiveGroups.has_count(64)
false
```
"""
has_count(deg::Int) = has(deg)

"""
    TransitiveGroups.has_identification(deg::Int)

Return whether identification of transitive groups of degree `deg` is available
via `TransitiveGroups.identification(G)`, where `G` is a permutation group
of degree `deg`.

# Examples
```jldoctest
julia> TransitiveGroups.has_identification(30)
true

julia> TransitiveGroups.has_identification(64)
false
```
"""
has_identification(deg::Int) = has(deg)

"""
    TransitiveGroups.has(deg::Int)

Return whether the transitive groups of degree `deg` are available for
use. This function should be used to test for the scope of the library
available.

# Examples
```jldoctest
julia> TransitiveGroups.has(30)
true

julia> TransitiveGroups.has(64)
false
```
"""
function has(deg::Int)
  @req deg >= 1 "degree must be positive, not $deg"
  return deg == 1 || GAP.Globals.TransitiveGroupsAvailable(deg)::Bool
end


"""
    TransitiveGroups.count(deg::Int)

Return the number of transitive groups of degree `deg`,
up to permutation isomorphism.

# Examples
```jldoctest
julia> TransitiveGroups.count(30)
5712

julia> TransitiveGroups.count(64)
ERROR: ArgumentError: the number of transitive groups of degree 64 is not available
[...]
```
"""
function count(deg::Int)
  @req has_count(deg) "the number of transitive groups of degree $(deg) is not available"
  return GAP.Globals.NrTransitiveGroups(deg)::Int
end

"""
    TransitiveGroups.get(deg::Int, i::Int)

Return the `i`-th group in the catalogue of transitive groups over the set
`{1, ..., deg}` in GAP's Transitive Groups Library.
The output is a group of type `PermGroup`.

# Examples
```jldoctest
julia> TransitiveGroups.get(5, 4)
Alternating group of degree 5

julia> TransitiveGroups.get(5, 6)
ERROR: ArgumentError: there are only 5 transitive groups of degree 5, not 6
[...]
```
"""
function get(deg::Int, i::Int)
  @req has(deg) "transitive groups of degree $deg are not available"
  N = count(deg)
  @req i <= N "there are only $N transitive groups of degree $deg, not $i"
  return PermGroup(GAP.Globals.TransitiveGroup(deg,i), deg)
end

"""
    TransitiveGroups.identification(G::PermGroup)

Return a pair `(d, n)` such that `G` is permutation isomorphic with
`TransitiveGroups.get(d, n)`, where `G` acts transitively on `d` points.

If `G` is not transitive on its moved points, or if the transitive groups of
degree `d` are not available, an exception is thrown.

# Examples
```jldoctest
julia> G = symmetric_group(7);  m = TransitiveGroups.identification(G)
(7, 7)

julia> order(TransitiveGroups.get(m...)) == order(G)
true

julia> S = sub(G, [perm([1, 3, 4, 5, 2])])[1]
Permutation group of degree 7

julia> is_transitive(S)
false

julia> is_transitive(S, moved_points(S))
true

julia> m = TransitiveGroups.identification(S)
(4, 1)

julia> order(TransitiveGroups.get(m...)) == order(S)
true

julia> TransitiveGroups.identification(symmetric_group(64))
ERROR: ArgumentError: identification of transitive groups of degree 64 are not available
[...]

julia> S = sub(G, [perm([1, 3, 4, 5, 2, 7, 6])])[1];

julia> TransitiveGroups.identification(S)
ERROR: ArgumentError: group is not transitive on its moved points
[...]
```
"""
function identification(G::PermGroup)
  moved = moved_points(G)
  @req is_transitive(G, moved) "group is not transitive on its moved points"
  deg = length(moved)
  @req has(deg) "identification of transitive groups of degree $(deg) are not available"
  res = GAP.Globals.TransitiveIdentification(GapObj(G))::Int
  return deg, res
end

"""
    TransitiveGroups.get_all(L...)
    TransitiveGroups.get_one(L...)

`TransitiveGroups.get_all` returns the list of all transitive groups
(up to permutation isomorphism) satisfying the conditions described
by the arguments.
`TransitiveGroups.get_one` returns one transitive group satisfying
the conditions if such a group exists, and `nothing` otherwise.
These conditions may be of one of the following forms:

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

The type of the returned group(s) is `PermGroup`.

# Examples
```jldoctest
julia> TransitiveGroups.get_all(4)
5-element Vector{PermGroup}:
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 4
 Alternating group of degree 4
 Symmetric group of degree 4

julia> TransitiveGroups.get_one(4)
Permutation group of degree 4

julia> TransitiveGroups.get_all(degree => 3:5, is_abelian)
4-element Vector{PermGroup}:
 Alternating group of degree 3
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 5

julia> TransitiveGroups.get_one(degree => 3:5, is_abelian)
Alternating group of degree 3
```
"""
function get_all(L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (degree => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L; filter_attrs = _permgroup_filter_attrs)
   K = GAP.Globals.AllTransitiveGroups(gapargs...)::GapObj
   return [PermGroup(x) for x in K]
end
# TODO: turn this into an iterator, perhaps using PrimitiveGroupsIterator

function get_one(L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (degree => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L; filter_attrs = _permgroup_filter_attrs)
   K = GAP.Globals.OneTransitiveGroup(gapargs...)::GapObj
   K === GAP.Globals.fail && return nothing
   return PermGroup(K)
end
@doc (@doc get_all) get_one

end # module TransitiveGroups

export TransitiveGroups
