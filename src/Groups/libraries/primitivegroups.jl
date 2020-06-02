export
    all_primitive_groups,
    blocks,
    isprimitive,
    maximal_blocks,
    number_primitive_groups,
    primitive_group,
    representatives_minimal_blocks
    


###################################################################
# Primitive groups, block system
###################################################################

"""
    blocks(G::PermGroup, L::AbstractVector{Int})
Return a block system for the action of `G` over `L`, i.e. a minimal non-trivial partition of `L` preserved by the action of `G`. Here, `L` must be a subvector of [1..deg(`G`)] and it is considered only the action of `H` on `L`, where `H` is the subgroup of `G` that moves only points in `L`. If this action is not transitive, then an ERROR is returned. If `L` is not specified, then `L` is taken as the set of moved points by `G`.
"""
function blocks(G::PermGroup, L::AbstractVector{Int})
   @assert istransitive(G,L) "The group action is not transitive"
   l = GAP.gap_to_julia(GAP.Globals.Blocks(G.X, GAP.julia_to_gap(L)))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

blocks(G::PermGroup) = blocks(G,[i for i in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))])

"""
    maximal_blocks(G::PermGroup, L::AbstractVector{Int})
Return a maximal block system for the action of `G` over `L`, i.e. a maximal non-trivial partition of `L` preserved by the action of `G`. Here, `L` must be a subvector of [1..deg(`G`)] and it is considered only the action of `H` on `L`, where `H` is the subgroup of `G` that moves only points in `L`. If this action is not transitive, then an ERROR is returned. If `L` is not specified, then `L` is taken as the set of moved points by `G`.
"""
function maximal_blocks(G::PermGroup, L::AbstractVector{Int})
   @assert istransitive(G,L) "The group action is not transitive"
   l = GAP.gap_to_julia(GAP.Globals.MaximalBlocks(G.X, GAP.julia_to_gap(L)))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

maximal_blocks(G::PermGroup) = maximal_blocks(G,[i for i in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))])

"""
    representatives_minimal_blocks(G::PermGroup, L::AbstractVector{Int})
Return a list of block representatives for all minimal non-trivial block systems for the action of `G` over `L`. Here, `L` must be a subvector of [1..deg(`G`)] and it is considered only the action of `H` on `L`, where `H` is the subgroup of `G` that moves only points in `L`. If this action is not transitive, then an ERROR is returned. If `L` is not specified, then `L` is taken as the set of moved points by `G`.
"""
function representatives_minimal_blocks(G::PermGroup, L::AbstractVector{Int})
   @assert istransitive(G,L) "The group action is not transitive"
   l = GAP.gap_to_julia(GAP.Globals.RepresentativesMinimalBlocks(G.X, GAP.julia_to_gap(L)))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

representatives_minimal_blocks(G::PermGroup) = minimal_blocks(G,[i for i in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))])

"""
    allblocks(G::PermGroup)
Return a list of representatives of all block systems for the action of `G` on the set of moved points of `G`.
"""
function allblocks(G::PermGroup)
   l = GAP.gap_to_julia(GAP.Globals.AllBlocks(G.X))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end
"""

    isprimitive(G::PermGroup, L::AbstractVector{Int})
Return whether the action of the group `G` on `L` is primitive. If `L` is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
isprimitive(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(L))
isprimitive(G::PermGroup) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(1:G.deg))


"""
    number_primitive_groups(n::Int)
Return the number of primitive groups acting on a set of size `n`.
"""
number_primitive_groups(n::Int) = GAP.Globals.NrPrimitiveGroups(n)

"""
    primitive_group(deg::Int, i::Int)
Return the `i`-th group in the catalogue of primitive groups over the set {`1`,...,`deg`} in the GAP Small Groups Library. The output is a group of type ``PermGroup``.
"""
function primitive_group(deg::Int, n::Int)
   @assert n<= number_primitive_groups(deg) "There are less than $n primitive groups of degree $deg."
   return PermGroup(GAP.Globals.PrimitiveGroup(deg,n), deg)
end

"""
    all_primitive_groups(L...)
Return the list of all primitive groups (up to isomorphism) satisfying the conditions in `L`. Here, `L` is a vector whose arguments are organized as `L` = [ `func1`, `arg1`, `func2`, `arg2`, ... ], and the function returns all the groups `G` satisfying the conditions `func1`(`G`) = `arg1`, `func2`(`G`) = `arg2`, etc. An argument can be omitted if it corresponds to the boolean value ``true``.

# Example
```jldoctest
julia> all_primitive_groups(degree, 4, isabelian)
0-element Array{PermGroup,1}
```
returns the list of all abelian primitive groups acting on a set of order 4.

The type of the groups is ``PermGroup``.
"""
function all_primitive_groups(L...)
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

   K = GAP.Globals.CallFuncList(GAP.Globals.AllPrimitiveGroups,L1)
   return [PermGroup(K[i]) for i in 1:length(K)]          # GAP.julia_to_gap(K) does not work
end


