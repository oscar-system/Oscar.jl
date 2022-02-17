export
    all_transitive_groups,
    number_transitive_groups,
    transitive_group,
    transitive_identification
    

"""
    number_transitive_groups(n::Int)

Return the number of transitive groups acting on a set of size `n`,
up to permutation isomorphism.
"""
number_transitive_groups(n::Int) = GAP.Globals.NrTransitiveGroups(n)

"""
    transitive_group(deg::Int, i::Int)

Return the `i`-th group in the catalogue of transitive groups over the set
`{1, ..., deg}` in GAP's Transitive Groups Library.
The output is a group of type `PermGroup`.
"""
function transitive_group(deg::Int, n::Int)
   N = number_transitive_groups(deg)
   @assert n <= N "There are only $N transitive groups of degree $deg, up to permutation isomorphism."

   return PermGroup(GAP.Globals.TransitiveGroup(deg,n), deg)
end

"""
    transitive_identification(G::PermGroup)

Return `m` such that `G` is permutation isomorphic with
`transitive_group(d, m)`, where `G` is transitive on `d` points,
with `d < 32`.

If `G` moves more than 31 points then `-1` is returned.
Otherwise, if `G` is not transitive on its moved points then `0` is returned,

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
-1

julia> S = sub(G, [gap_perm([1,3,4,5,2,7,6])])[1];

julia> istransitive(S, moved_points(S))
false

julia> transitive_identification(S)
0
```
"""
function transitive_identification(G::PermGroup)
  moved = moved_points(G)
  length(moved) < 32 || return -1
  istransitive(G, moved) || return 0
  return GAP.Globals.TransitiveIdentification(G.X)
end


"""
    all_transitive_groups(L...)

Return the list of all transitive groups (up to permutation isomorphism)
satisfying the conditions in `L`.
Here, `L` is a vector whose arguments are organized as `L` = [ `func1`,
`arg1`, `func2`, `arg2`, ... ], and the function returns all the groups `G`
satisfying the conditions `func1`(`G`) = `arg1`, `func2`(`G`) = `arg2`, etc.
An argument can be omitted if it corresponds to the boolean value `true`.

# Examples
```jldoctest
julia> all_transitive_groups(degree, 4, isabelian)
2-element Vector{PermGroup}:
 C(4) = 4
 E(4) = 2[x]2
```
returns the list of all abelian transitive groups acting on a set of order 4.

The type of the groups is `PermGroup`.
"""
function all_transitive_groups(L...)
   valid, temp = CheckValidType(L; isapg=true)
   @assert valid "Wrong type inserted"
   isargument = false                     # says if the inserted value is the argument of the previous value
   
   L1 = Vector(undef, length(L)+temp)
   pos = 1
   for i in 1:length(L)
      if typeof(L[i]) <: Function
         if isargument
            L1[pos] = true
            pos += 1
         end
         L1[pos] = find_index_function(L[i], true)[2]
         isargument = true
      else
         L1[pos] = GAP.julia_to_gap(L[i])
         isargument = false
      end
   pos+=1
   end
   if isargument L1[length(L1)]=true end

   K = GAP.Globals.AllTransitiveGroups(L1...)
   return [PermGroup(x) for x in K]          # GAP.julia_to_gap(K) does not work
end

