# This file contains code related to G-sets

if isdefined(Hecke, :isregular)
    import Hecke: isregular
end

export
    blocks,
    isprimitive,
    isregular,
    issemiregular,
    istransitive,
    maximal_blocks,
    representatives_minimal_blocks,
    transitivity

export gset, orbits, as_gset

# G-sets are "sets" (in a very general sense, these do not need to be objects of type `Set`)
# with an action by a group G::T.
# Alternatively one can define G-sets as a union of G-orbits.
# Potential examples include:
# - orbits of integers under a permutation
# - conjugacy classes of group elements
# - conjugacy classes of subgroups
# - block system
abstract type GSet{TG,TE} end


# TODO: add lots of concrete subtypes constructors, e.g. for
# - regular action of a group on itself
# - action of a perm group on its moved points
# - ...

# G-set that is willing to write down complete orbits and elements set
mutable struct GSet_By_Elements{TG,TE} <: GSet{TG,TE}
    group::TG
    action_function::Function
    seeds::Set{TE}
    AbstractAlgebra.@declare_other  # we want to store the orbits, the set of elements etc. if known

    function GSet_By_Elements(G::TG, fun::Function, seeds::Set{TE}) where TG <: GAPGroup where TE
        @assert length(seeds) > 0
        return new{TG,TE}(G, fun, seeds, Dict{Symbol,Any}())
    end
end

#############################################################################
##
##  general method with explicit action function

gset(G::GAPGroup, fun::Function, Omega::Set) = GSet_By_Elements(G, fun, Omega)


#############################################################################
##
##  G-sets where the action function can be omitted

## natural action of permutations on positive integers
gset(G::PermGroup, Omega::Set{T}) where T<:Union{Base.Integer,fmpz} = GSet_By_Elements(G, ^, Omega)

## action of permutations on sets of positive integers
gset(G::PermGroup, Omega::Set{Set{T}}) where T<:Union{Base.Integer,fmpz} = GSet_By_Elements(G, on_sets, Omega)

## action of permutations on tuples of positive integers
## (given either as vectors or as tuples)
gset(G::PermGroup, Omega::Set{Vector{T}}) where T<:Union{Base.Integer,fmpz} = GSet_By_Elements(G, on_tuples, Omega)

gset(G::PermGroup, Omega::Set{Tuple{Vararg{T}}}) where T<:Union{Base.Integer,fmpz} = GSet_By_Elements(G, on_tuples, Omega)

function gset(G::PermGroup)
    omega = gset(G, Set(1:G.deg))
    AbstractAlgebra.set_special(omega, :elements => omega.seeds)
    return omega
end


#############################################################################
##
##  G-sets given by the complete set

function as_gset(G::TG, fun::Function, Omega::Set) where TG <: GAPGroup
    omega = GSet_By_Elements(G, fun, Omega)
    AbstractAlgebra.set_special(omega, :elements => omega.seeds)
    return omega
end

as_gset(G::TG, Omega::Set) where TG <: GAPGroup = as_gset(G, ^, Omega)


#############################################################################
##
##  wrapper objects for elements of G-sets,
##  with fields `gset` (the G-set) and `objects` (the unwrapped object),
##  are used to
##  - apply group elements via `^`,
##    not via the action function stored in the G-set
##  - write something like `stabilizer(omega)`

struct ElementOfGSet
    gset::GSet
    obj::Any
end

function (Omega::GSet)(obj::Any)
    return ElementOfGSet(Omega, obj)
end

function ^(omega::ElementOfGSet, g::T) where {T<:AbstractAlgebra.GroupElem}
    Omega = omega.gset
    return ElementOfGSet(Omega, Omega.action_function(omega.obj, g))
end

==(omega1::ElementOfGSet, omega2::ElementOfGSet) = ((omega1.gset == omega2.gset) && (omega1.obj == omega2.obj))

orbit(omega::ElementOfGSet) = orbit(omega.gset, omega.obj)


#############################################################################
##
##  `:orbit` a Set of points

function orbit(Omega::GSet_By_Elements{<:GAPGroup,T}, omega) where T
    acts = gens(Omega.group)
    orbarray = [omega]
    orb = Set(orbarray)
    fun = Omega.action_function
    for p in orbarray
      for g in acts
        img = fun(p, g)
        if !(img in orb)
          push!(orbarray, img)
          push!(orb, img)
        end
      end
    end
    return orb
end

orbit(G::PermGroup, pt::T) where {T<:Union{Base.Integer,fmpz}} = GSet_By_Elements(G, Set(pt))


#############################################################################
##
##  `:orbits` a Set of Sets of points

function orbits(Omega::GSet_By_Elements{<:GAPGroup,T}) where T
    orbs = AbstractAlgebra.get_special(Omega, :orbits)
    if orbs == nothing
      G = Omega.group
      orbs = Set{typeof(Omega.seeds)}()
      for p in Omega.seeds
        if all(o -> !(p in o), orbs)
          push!(orbs, orbit(Omega, p))
        end
      end
      AbstractAlgebra.set_special(Omega, :orbits => orbs)
    end
    return orbs
end


#############################################################################
##
##  `:elements` a Set of points

function elements(Omega::GSet_By_Elements)
    elms = AbstractAlgebra.get_special(Omega, :elements)
    if elms == nothing
      orbs = orbits(Omega)
      elms = union(orbs...)
      AbstractAlgebra.set_special(Omega, :elements => elms)
    end
    return elms
end

Base.length(Omega::GSet) = length(elements(Omega))

representative(Omega::GSet) = first(Omega.seeds)

acting_domain(Omega::GSet) = Omega.group

# TODO:
# - Base.iterate


blocks(G::GSet) = error("not implemented")
maximal_blocks(G::GSet) = error("not implemented")
representatives_minimal_blocks(G::GSet) = error("not implemented")
allblocks(G::GSet) = error("not implemented")

function istransitive(Omega::GSet)
    length(Omega.seeds) == 1 && return true
    return length(orbits(Omega)) == 1
end

isregular(Omega::GSet) = istransitive(Omega) && length(Omega) == order(Omega.group)

function issemiregular(Omega::GSet)
    ord = order(Omega.group)
    return all(orb -> length(orb) == ord, orbits(Omega))
end

isprimitive(G::GSet) = error("not implemented")


"""
    blocks(G::PermGroup, L::AbstractVector{Int})

Return a block system for the action of `G` over `L`, i.e. a minimal
non-trivial partition of `L` preserved by the action of `G`. Here, `L`
must be a subvector of [1..degree(`G`)] and it is considered only the
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
Here, `L` must be a subvector of [1..degree(`G`)] and it is considered only
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
[1..degree(`G`)] and it is considered only the action of `H` on `L`, where
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
not specified, then `L` is taken as [1,...,degree(`G`)].
"""
transitivity(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.Transitivity(G.X, GAP.julia_to_gap(L))
transitivity(G::PermGroup) = GAP.Globals.Transitivity(G.X, GAP.julia_to_gap(1:G.deg))

"""
    istransitive(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is transitive. If `L`
is not specified, then `L` is taken as [1,...,degree(`G`)].
"""
istransitive(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(L))
istransitive(G::PermGroup) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(1:G.deg))

"""
    isprimitive(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is primitive. If `L`
is not specified, then `L` is taken as [1,...,degree(`G`)].
"""
isprimitive(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(L))
isprimitive(G::PermGroup) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(1:G.deg))

"""
    isregular(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is regular (i.e.
transitive and semiregular). If `L` is not specified, then `L` is taken
as [1,...,degree(`G`)].
"""
isregular(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsRegular(G.X, GAP.julia_to_gap(L))
isregular(G::PermGroup) = GAP.Globals.IsRegular(G.X, GAP.julia_to_gap(1:G.deg))

"""
    issemiregular(G::PermGroup, L::AbstractVector{Int})

Return whether the action of the group `G` on `L` is semiregular (i.e.
the stabilizer of each point is the identity). If `L` is not specified,
then `L` is taken as [1,...,degree(`G`)].
"""
issemiregular(G::PermGroup, L::AbstractVector{Int}) = GAP.Globals.IsSemiRegular(G.X, GAP.julia_to_gap(L))
issemiregular(G::PermGroup) = GAP.Globals.IsSemiRegular(G.X, GAP.julia_to_gap(1:G.deg))

