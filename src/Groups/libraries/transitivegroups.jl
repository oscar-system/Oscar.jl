export
    all_transitive_groups,
    istransitive,
    isregular,
    issemiregular,
    number_transitive_groups,
    transitive_group,
    transitive_identification,
    transitivity
    

"""
    number_transitive_groups(n::Int)
Return the number of transitive groups acting on a set of size `n`.
"""
number_transitive_groups(n::Int) = GAP.Globals.NrTransitiveGroups(n)

"""
    transitive_group(deg::Int, i::Int)
Return the `i`-th group in the catalogue of transitive groups over the set {`1`,...,`deg`} in the GAP Small Groups Library. The output is a group of type ``PermGroup``.
"""
function transitive_group(deg::Int, n::Int)
   @assert n<= number_transitive_groups(deg) "There are less than $n transitive groups of degree $deg."
   return PermGroup(GAP.Globals.TransitiveGroup(deg,n), deg)
end

"""
    transitive_identification(G::PermGroup)
Return (`deg`, `m`), where `G` = transitive_group(`deg`,`m`).
"""
transitive_identification(G::PermGroup) = GAP.Globals.TransitiveIdentification(G.X)

"""
    transitivity(G::PermGroup, L::AbstractVector{Int})
Return the maximum `k` such that the action of `G` over `L` is `k`-transitive. The output is ``0`` if `G` is not transitive. If `L` is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
transitivity(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.Transitivity(G.X, GAP.julia_to_gap(L))
transitivity(G::PermGroup) = GAP.Globals.Transitivity(G.X, GAP.julia_to_gap(1:G.deg))

"""
    istransitive(G::PermGroup, L::AbstractVector{Int})
Return whether the action of the group `G` on `L` is transitive. If `L` is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
istransitive(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(L))
istransitive(G::PermGroup) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(1:G.deg))

"""
    isregular(G::PermGroup, L::AbstractVector{Int})
Return whether the action of the group `G` on `L` is regular (i.e. transitive and semiregular). If `L` is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
isregular(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsRegular(G.X, GAP.julia_to_gap(L))
isregular(G::PermGroup) = GAP.Globals.IsRegular(G.X, GAP.julia_to_gap(1:G.deg))

"""
    issemiregular(G::PermGroup, L::AbstractVector{Int})
Return whether the action of the group `G` on `L` is semiregular (i.e. the stabilizer of each point is the identity). If `L` is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
issemiregular(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsSemiRegular(G.X, GAP.julia_to_gap(L))
issemiregular(G::PermGroup) = GAP.Globals.IsSemiRegular(G.X, GAP.julia_to_gap(1:G.deg))


"""
    all_transitive_groups(L...)
Return the list of all transitive groups (up to isomorphism) satisfying the conditions in `L`. Here, `L` is a vector whose arguments are organized as `L` = [ `func1`, `arg1`, `func2`, `arg2`, ... ], and the function returns all the groups `G` satisfying the conditions `func1`(`G`) = `arg1`, `func2`(`G`) = `arg2`, etc. An argument can be omitted if it corresponds to the boolean value ``true``.

# Example
```jldoctest
julia> all_transitive_groups(degree, 4, isabelian)
2-element Array{PermGroup,1}:
 C(4) = 4
 E(4) = 2[x]2
```
returns the list of all abelian transitive groups acting on a set of order 4.

The type of the groups is ``PermGroup``.
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
   L1 = GAP.julia_to_gap(L1)

   K = GAP.Globals.CallFuncList(GAP.Globals.AllTransitiveGroups,L1)
   return [PermGroup(K[i]) for i in 1:length(K)]          # GAP.julia_to_gap(K) does not work
end

