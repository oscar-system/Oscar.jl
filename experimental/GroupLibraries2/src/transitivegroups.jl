"""
    TransitiveGroupsLibrary

This type is intended to be used as the first argument of functions
such as `get`, `count`, `identification` in the context of
GAP's library of transitive permutation groups.

This library provides access to the transitive permutation groups
up to degree 48, via the GAP package `TransGrp` [Hul23](@cite).

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
struct TransitiveGroupsLibrary end

"""
    has_count(TransitiveGroupsLibrary, deg::Int)

Return whether the number of transitive groups of degree `deg` is available for
use via `count(TransitiveGroupsLibrary, deg)`.

# Examples
```jldoctest
julia> has_count(TransitiveGroupsLibrary, 30)
true

julia> has_count(TransitiveGroupsLibrary, 64)
false
```
"""
has_count(T::Type{TransitiveGroupsLibrary}, deg::Int) = has(T, deg)

"""
    has_identification(TransitiveGroupsLibrary, deg::Int)

Return whether identification of transitive groups of degree `deg` is available
via `identification(TransitiveGroupsLibrary, G)`, where `G` is a permutation group
of degree `deg`.

# Examples
```jldoctest
julia> has_identification(TransitiveGroupsLibrary, 30)
true

julia> has_identification(TransitiveGroupsLibrary, 64)
false
```
"""
has_identification(T::Type{TransitiveGroupsLibrary}, deg::Int) = has(T, deg)

"""
    has(TransitiveGroupsLibrary, deg::Int)

Return whether the transitive groups of degree `deg` are available for
use. This function should be used to test for the scope of the library
available.

# Examples
```jldoctest
julia> has(TransitiveGroupsLibrary, 30)
true

julia> has(TransitiveGroupsLibrary, 64)
false
```
"""
function has(::Type{TransitiveGroupsLibrary}, deg::Int)
  @req deg >= 1 "degree must be positive, not $deg"
  return deg == 1 || GAP.Globals.TransitiveGroupsAvailable(deg)::Bool
end


"""
    count(TransitiveGroupsLibrary, deg::Int)

Return the number of transitive groups of degree `deg`,
up to permutation isomorphism.

# Examples
```jldoctest
julia> count(TransitiveGroupsLibrary, 30)
5712

julia> count(TransitiveGroupsLibrary, 64)
ERROR: ArgumentError: the number of transitive groups of degree 64 is not available
[...]
```
"""
function Base.count(T::Type{TransitiveGroupsLibrary}, deg::Int)
  @req has_count(T, deg) "the number of transitive groups of degree $(deg) is not available"
  return GAP.Globals.NrTransitiveGroups(deg)::Int
end

"""
    get(TransitiveGroupsLibrary, deg::Int, i::Int)

Return the `i`-th group in the catalogue of transitive groups over the set
`{1, ..., deg}` in GAP's Transitive Groups Library.
The output is a group of type `PermGroup`.

# Examples
```jldoctest
julia> get(TransitiveGroupsLibrary, 5, 4)
Alternating group of degree 5

julia> get(TransitiveGroupsLibrary, 5, 6)
ERROR: ArgumentError: there are only 5 transitive groups of degree 5, not 6
[...]
```
"""
function Base.get(T::Type{TransitiveGroupsLibrary}, deg::Int, i::Int)
  @req has(T, deg) "transitive groups of degree $deg are not available"
  N = count(T, deg)
  @req i <= N "there are only $N transitive groups of degree $deg, not $i"
  return PermGroup(GAP.Globals.TransitiveGroup(deg,i), deg)
end

"""
    identification(TransitiveGroupsLibrary, G::PermGroup)

Return a pair `(d, n)` such that `G` is permutation isomorphic with
`get(TransitiveGroupsLibrary, d, n)`, where `G` acts transitively on `d` points.

If `G` is not transitive on its moved points, or if the transitive groups of
degree `d` are not available, an exception is thrown.

# Examples
```jldoctest
julia> G = symmetric_group(7);  m = identification(TransitiveGroupsLibrary, G)
(7, 7)

julia> order(get(TransitiveGroupsLibrary, m...)) == order(G)
true

julia> S = sub(G, [perm([1, 3, 4, 5, 2])])[1]
Permutation group of degree 7

julia> is_transitive(S)
false

julia> is_transitive(S, moved_points(S))
true

julia> m = identification(TransitiveGroupsLibrary, S)
(4, 1)

julia> order(get(TransitiveGroupsLibrary, m...)) == order(S)
true

julia> identification(TransitiveGroupsLibrary, symmetric_group(64))
ERROR: ArgumentError: identification of transitive groups of degree 64 are not available
[...]

julia> S = sub(G, [perm([1, 3, 4, 5, 2, 7, 6])])[1];

julia> identification(TransitiveGroupsLibrary, S)
ERROR: ArgumentError: group is not transitive on its moved points
[...]
```
"""
function identification(T::Type{TransitiveGroupsLibrary}, G::PermGroup)
  moved = moved_points(G)
  @req is_transitive(G, moved) "group is not transitive on its moved points"
  deg = length(moved)
  @req has(T, deg) "identification of transitive groups of degree $(deg) are not available"
  res = GAP.Globals.TransitiveIdentification(GapObj(G))::Int
  return deg, res
end

_docstr = """
    get_all(TransitiveGroupsLibrary, L...)
    get_one(TransitiveGroupsLibrary, L...)

`get_all` returns the list of all transitive groups
(up to permutation isomorphism) satisfying the conditions described
by the arguments.

`get_one` returns one transitive group satisfying
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
julia> get_all(TransitiveGroupsLibrary, 4)
5-element Vector{PermGroup}:
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 4
 Alternating group of degree 4
 Symmetric group of degree 4

julia> get_one(TransitiveGroupsLibrary, 4)
Permutation group of degree 4

julia> get_all(TransitiveGroupsLibrary, degree => 3:5, is_abelian)
4-element Vector{PermGroup}:
 Alternating group of degree 3
 Permutation group of degree 4
 Permutation group of degree 4
 Permutation group of degree 5

julia> get_one(TransitiveGroupsLibrary, degree => 3:5, is_abelian)
Alternating group of degree 3
```
"""

@doc _docstr
function get_all(::Type{TransitiveGroupsLibrary}, L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (degree => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L; filter_attrs = _permgroup_filter_attrs)
   K = GAP.Globals.AllTransitiveGroups(gapargs...)::GapObj
   return [PermGroup(x) for x in K]
end

@doc _docstr
function get_one(::Type{TransitiveGroupsLibrary}, L...)
   @req !isempty(L) "must specify at least one filter"
   if L[1] isa IntegerUnion || L[1] isa AbstractVector{<:IntegerUnion}
      L = (degree => L[1], L[2:end]...)
   end
   gapargs = translate_group_library_args(L; filter_attrs = _permgroup_filter_attrs)
   K = GAP.Globals.OneTransitiveGroup(gapargs...)::GapObj
   K === GAP.Globals.fail && return nothing
   return PermGroup(K)
end
