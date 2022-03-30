export
    all_small_groups,
    number_small_groups,
    small_group,
    small_group_identification
    


###################################################################
# Small groups
###################################################################

"""
    small_group(n::Int, i::Int)

Return the `i`-th group of order `n` in the catalogue of GAP's
Small Groups Library.
The group is given of type `PcGroup` if the group is solvable,
`PermGroup` otherwise.
"""
function small_group(n::Int, m::Int)
  N = number_small_groups(n)
  @assert m <= N "There are only $N groups of order $n, up to isomorphism."
  G = GAP.Globals.SmallGroup(n, m)
  T = _get_type(G)
  return T(G)
end

"""
    small_group_identification(G::Group)

Return `(n, m)`, where `G` is isomorphic with `small_group(n, m)`.
"""
function small_group_identification(G::GAPGroup)
  r = GAP.Globals.IdGroup(G.X)
  return (r[1], r[2])
end

"""
    number_small_groups(n::Int)

Return the number of groups of order `n`, up to isomorphism.
"""
number_small_groups(n::Int) = GAP.Globals.NumberSmallGroups(n)


"""
    all_small_groups(n::Int, L...)

Return the list of all groups (up to isomorphism) of order `n` and satisfying
the conditions in `L`. Here, `L` is a vector whose arguments are organized as
`L` = [ `func1`, `arg1`, `func2`, `arg2`, ... ], and the function returns all
the groups `G` satisfying the conditions `func1`(`G`) = `arg1`, `func2`(`G`) =
`arg2`, etc. An argument can be omitted if it corresponds to the boolean value
`true`.

The following command returns the list of all abelian non-cyclic groups
of order 12.

# Examples
```jldoctest
julia> all_small_groups(12, iscyclic, false, isabelian)
1-element Vector{PcGroup}:
 <pc group of size 12 with 3 generators>

```

The type of the groups is `PcGroup` if the group is solvable, `PermGroup` otherwise.
"""
function all_small_groups(n::Int, L...)
   @assert CheckValidType(L)[1] "Wrong type inserted"
   
   L1 = Vector(undef, length(L)+1)
   L1[1] = n
   for i in 1:length(L)
      if typeof(L[i]) <: Function
         L1[i+1] = find_index_function(L[i],false)[2]
      else
         L1[i+1] = GAP.julia_to_gap(L[i])
      end
   end

   K = GAP.Globals.AllSmallGroups(L1...)
   return [_get_type(x)(x) for x in K]          # GAP.julia_to_gap(K) does not work
end
#T what does this comment mean?

#T problem:

#T all_small_groups( 60, issimple ) -> Array{PermGroup,1}
#T all_small_groups( 60 )  -> Array{Oscar.GAPGroup,1}
#T all_small_groups( 59 )  -> Array{PcGroup,1}

#T Do we want this???
