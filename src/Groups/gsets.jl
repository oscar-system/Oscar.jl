# This file contains code related to G-sets

export
    blocks,
    isprimitive,
    isregular,
    issemiregular,
    istransitive,
    maximal_blocks,
    representatives_minimal_blocks,
    transitivity


# G-sets are "sets" (in a very general sense, these do not need to be objects of type `Set`)
# with an action by a group G::T.
# Alternatively one can define G-sets as a union of G-orbits.
# Potential examples include:
# - orbits of integers under a permutation
# - conjugacy classes of group elements
# - conjugacy classes of subgroups
# - block system
abstract type GSet{T} end

# TODO: add lots of concrete subtypes constructors, e.g. for
# - regular action of a group on itself
# - action of a perm group on its moved points
# - ...


Base.length(C::GSet) = error("not implemented")
representative(C::GSet) = error("not implemented")
elements(C::GSet) = error("not implemented")
acting_domain(C::GSet) = error("not implemented")
# TODO:
# - Base.iterate


blocks(G::GSet) = error("not implemented")
maximal_blocks(G::GSet) = error("not implemented")
representatives_minimal_blocks(G::GSet) = error("not implemented")
allblocks(G::GSet) = error("not implemented")

istransitive(G::GSet) = error("not implemented")
isprimitive(G::GSet) = error("not implemented")
isregular(G::GSet) = error("not implemented")
issemiregular(G::GSet) = error("not implemented")


"""
    blocks(G::PermGroup, L::AbstractVector{Int})

Return a block system for the action of `G` over `L`, i.e. a minimal
non-trivial partition of `L` preserved by the action of `G`. Here, `L`
must be a subvector of [1..deg(`G`)] and it is considered only the
action of `H` on `L`, where `H` is the subgroup of `G` that moves only
points in `L`. If this action is not transitive, then an ERROR is
returned. If `L` is not specified, then `L` is taken as the set of moved
points by `G`.
"""
function blocks(G::PermGroup, L::AbstractVector{Int})
   @assert istransitive(G,L) "The group action is not transitive"
   l = GAP.gap_to_julia(GAP.Globals.Blocks(G.X, GAP.julia_to_gap(L)))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

blocks(G::PermGroup) = blocks(G,[i for i in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))])

"""
    maximal_blocks(G::PermGroup, L::AbstractVector{Int})

Return a maximal block system for the action of `G` over `L`, i.e. a
maximal non-trivial partition of `L` preserved by the action of `G`.
Here, `L` must be a subvector of [1..deg(`G`)] and it is considered only
the action of `H` on `L`, where `H` is the subgroup of `G` that moves
only points in `L`. If this action is not transitive, then an ERROR is
returned. If `L` is not specified, then `L` is taken as the set of moved
points by `G`.
"""
function maximal_blocks(G::PermGroup, L::AbstractVector{Int})
   @assert istransitive(G,L) "The group action is not transitive"
   l = GAP.gap_to_julia(GAP.Globals.MaximalBlocks(G.X, GAP.julia_to_gap(L)))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

maximal_blocks(G::PermGroup) = maximal_blocks(G,[i for i in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))])

"""
    representatives_minimal_blocks(G::PermGroup, L::AbstractVector{Int})

Return a list of block representatives for all minimal non-trivial block
systems for the action of `G` over `L`. Here, `L` must be a subvector of
[1..deg(`G`)] and it is considered only the action of `H` on `L`, where
`H` is the subgroup of `G` that moves only points in `L`. If this action
is not transitive, then an ERROR is returned. If `L` is not specified,
then `L` is taken as the set of moved points by `G`.
"""
function representatives_minimal_blocks(G::PermGroup, L::AbstractVector{Int})
   @assert istransitive(G,L) "The group action is not transitive"
   l = GAP.gap_to_julia(GAP.Globals.RepresentativesMinimalBlocks(G.X, GAP.julia_to_gap(L)))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

representatives_minimal_blocks(G::PermGroup) = minimal_blocks(G,[i for i in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))])

"""
    allblocks(G::PermGroup)

Return a list of representatives of all block systems for the action of
`G` on the set of moved points of `G`.
"""
function allblocks(G::PermGroup)
   l = GAP.gap_to_julia(GAP.Globals.AllBlocks(G.X))
   return [ [y for y in l1] for l1 in l]              # to return an array of integers
end

"""
    transitivity(G::PermGroup, L::AbstractVector{Int})

Return the maximum `k` such that the action of `G` over `L` is
`k`-transitive. The output is ``0`` if `G` is not transitive. If `L` is
not specified, then `L` is taken as [1,...,deg(`G`)].
"""
transitivity(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.Transitivity(G.X, GAP.julia_to_gap(L))
transitivity(G::PermGroup) = GAP.Globals.Transitivity(G.X, GAP.julia_to_gap(1:G.deg))

"""
    istransitive(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is transitive. If `L`
is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
istransitive(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(L))
istransitive(G::PermGroup) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(1:G.deg))

"""
    isprimitive(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is primitive. If `L`
is not specified, then `L` is taken as [1,...,deg(`G`)].
"""
isprimitive(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(L))
isprimitive(G::PermGroup) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(1:G.deg))

"""
    isregular(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is regular (i.e.
transitive and semiregular). If `L` is not specified, then `L` is taken
as [1,...,deg(`G`)].
"""
isregular(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsRegular(G.X, GAP.julia_to_gap(L))
isregular(G::PermGroup) = GAP.Globals.IsRegular(G.X, GAP.julia_to_gap(1:G.deg))

"""
    issemiregular(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is semiregular (i.e.
the stabilizer of each point is the identity). If `L` is not specified,
then `L` is taken as [1,...,deg(`G`)].
"""
issemiregular(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsSemiRegular(G.X, GAP.julia_to_gap(L))
issemiregular(G::PermGroup) = GAP.Globals.IsSemiRegular(G.X, GAP.julia_to_gap(1:G.deg))

