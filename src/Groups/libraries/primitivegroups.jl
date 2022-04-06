export
    all_primitive_groups,
    has_primitive_groups,
    number_primitive_groups,
    primitive_group

# TODO: use PrimitiveGroupsIterator

###################################################################
# Primitive groups, block system
###################################################################

"""
    has_primitive_groups(n::Int)

Return `true` if the primitive permutation groups of degree `n` are available
via `primitive_group` and `all_primitive_groups`, otherwise `false`.

# Examples
```jldoctest
julia> has_primitive_groups(50)
true

julia> has_primitive_groups(5000)
false
```
"""
function has_primitive_groups(n::Int)
    n >= 1 || throw(ArgumentError("degree must be positive, not $n"))
    return GAP.Globals.PrimitiveGroupsAvailable(GAP.Obj(n))
end

"""
    number_primitive_groups(n::Int)

Return the number of primitive permutation groups of degree `n`,
up to permutation isomorphism.
"""
function number_primitive_groups(n::Int)
    has_primitive_groups(n) || throw(ArgumentError("primitive groups of degree $n are not currently available"))
    return fmpz(GAP.Globals.NrPrimitiveGroups(n))
end

"""
    primitive_group(deg::Int, i::Int)

Return the `i`-th group in the catalogue of primitive groups over the set
{`1`,...,`deg`} in the GAP Small Groups Library. The output is a group of type
``PermGroup``.
"""
function primitive_group(deg::Int, n::Int)
   has_primitive_groups(deg) || throw(ArgumentError("primitive groups of degree $deg are not currently available"))
   N = number_primitive_groups(deg)
   @assert n <= N "There are only $N primitive groups of degree $deg."
   return PermGroup(GAP.Globals.PrimitiveGroup(deg,n), deg)
end

"""
    all_primitive_groups(L...)

Return the list of all primitive groups (up to permutation isomorphism)
satisfying the conditions describe by the arguments. These conditions
may be of one of the following forms:

- `func => intval` selects groups for which the function `func` returns `intval`
- `func => list` selects groups for which the function `func` returns any element inside `list`
- `predicate` selects groups for which the function `predicate` returns `true`
- `!predicate` selects groups for which the function `predicate` returns `false`

The type of the returned groups is `PermGroup`.

# Examples
```jldoctest
julia> all_primitive_groups(degree => 3:5, isabelian)
2-element Vector{PermGroup}:
 A(3)
 C(5)
```
returns the list of all abelian primitive groups acting on 3, 4 or 5 points.
"""
function all_primitive_groups(L...)
   !isempty(L) || throw(ArgumentError("must specify at least one filter"))
   gapargs = translate_group_library_args(L; permgroups=true)
   K = GAP.Globals.AllPrimitiveGroups(gapargs...)
   return [PermGroup(x) for x in K]
end
