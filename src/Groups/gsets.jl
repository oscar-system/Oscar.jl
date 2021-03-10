# This file contains code related to G-sets
# The idea of the implementation is that available GAP functionality
# for computing orbits, permutations, actions, etc. can be used,
# but that local replacements by pure Julia code (independent of GAP)
# are welcome.

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

export GSet, gset, orbits, as_gset, unwrap, permutation, action_homomorphism,
    representative_action

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


#############################################################################
##
##  GSetByElements:
##  a G-set that is willing to write down complete orbits and elements lists;
##  fields are
##  - the group that acts, of type `T`,
##  - the Julia function (like `on_tuples`) that describes the action,
##  - the seeds (something iterable) whose closure under the action is the G-set
##  - the dictionary used to store attributes (orbits, elements, ...)

mutable struct GSetByElements{T} <: GSet{T}
    group::T
    action_function::Function
    seeds
    AbstractAlgebra.@declare_other

    function GSetByElements(G::T, fun::Function, seeds) where T<:GAPGroup
        @assert ! isempty(seeds)
        return new{T}(G, fun, seeds, Dict{Symbol,Any}())
    end
end
#TODO: How can I specify that `seeds` should be an iterable object?


#############################################################################
##
##  general method with explicit action function

gset(G::GAPGroup, fun::Function, Omega) = GSetByElements(G, fun, Omega)


#############################################################################
##
##  G-sets where the action function can be omitted
##
##  (We use an indirection via `gset_by_type`, in order to admit specifying
##  a default action depending on the element type of `Omega` (which can be
##  any iterable collection.)

gset(G::T, Omega) where T<:GAPGroup = gset_by_type(G, Omega, eltype(Omega))


## natural action of permutations on positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Union{Base.Integer,fmpz} = GSetByElements(G, ^, Omega)

## action of permutations on sets of positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Set{T2} where T2<:Union{Base.Integer,fmpz} = GSetByElements(G, on_sets, Omega)

## action of permutations on vectors of positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Vector{T2} where T2<:Union{Base.Integer,fmpz} = GSetByElements(G, on_tuples, Omega)

## action of permutations on tuples of positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Tuple{Vararg{T2}} where T2<:Union{Base.Integer,fmpz} = GSetByElements(G, on_tuples, Omega)

## (add more such actions: on sets of sets, on sets of tuples, ...)

## natural action of a permutation group on the integers 1, ..., degree
function gset(G::PermGroup)
    omega = gset(G, 1:G.deg)
    AbstractAlgebra.set_special(omega, :elements => omega.seeds)
    return omega
end


#############################################################################
##
#TODO: Compute membership without writing down all elements,
#      using what is called `RepresentativeAction` in GAP.

function Base.in(omega, Omega::GSet)
    omega in Omega.seeds && return true
    return omega in elements(Omega)
end


#############################################################################
##
##  G-sets given by the complete set

function as_gset(G::T, fun::Function, Omega) where T<:GAPGroup
    omega = GSetByElements(G, fun, Omega)
    AbstractAlgebra.set_special(omega, :elements => omega.seeds)
    return omega
end

as_gset(G::T, Omega) where T<:GAPGroup = as_gset(G, ^, Omega)


#############################################################################
##
##  wrapper objects for elements of G-sets,
##  with fields `gset` (the G-set) and `objects` (the unwrapped object)
##
##  These objects are optional ("syntactic sugar"), they can be used to
##  - apply group elements via `^`,
##    not via the action function stored in the G-set,
##  - write something like `orbit(omega)`, `stabilizer(omega)`.

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

Base.in(omega::ElementOfGSet, Omega::GSet) = Base.in(omega.obj, Omega)

Hecke.orbit(omega::ElementOfGSet) = Hecke.orbit(omega.gset, omega.obj)


unwrap(omega::Any) = omega

unwrap(omega::ElementOfGSet) = omega.obj


#############################################################################
##
##  `:orbit`

function Hecke.orbit(Omega::GSetByElements{<:GAPGroup}, omega)
    G = Omega.group
    acts = GAP.julia_to_gap(gens(G))
    gfun = GAP.julia_to_gap(Omega.action_function)

    # The following works only because GAP does not check
    # whether the given (dummy) group 'G.X' fits to the given generators,
    # or whether the elements of 'acts' are group elements.
    orb = GAP.gap_to_julia(GAP.Globals.Orbit(G.X, omega, acts, acts, gfun))

    res = as_gset(Omega.group, Omega.action_function, orb)
    # We know that this G-set is transitive.
    AbstractAlgebra.set_special(res, :orbits => [orb])
    return res
end

# simpleminded alternative directly in Julia
# In fact, '<:GAPGroup' is not used at all in this function.
function orbit_via_Julia(Omega::GSetByElements{<:GAPGroup}, omega)
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

    res = as_gset(Omega.group, Omega.action_function, orbarray)
    # We know that this G-set is transitive.
    AbstractAlgebra.set_special(res, :orbits => [orbarray])
    return res
end


#TODO: - Currently this conflicts with a method in
#        `experimental/GaloisGrp/Group.jl`
#      - Provide more convenience methods for "natural" actions",
#        as for `gset`?
# Hecke.orbit(G::PermGroup, pt::T) where {T<:Union{Base.Integer,fmpz}} = GSetByElements(G, [pt])


#############################################################################
##
##  `:orbits` an array of G-sets

function orbits(Omega::GSetByElements{<:GAPGroup})
    orbs = AbstractAlgebra.get_special(Omega, :orbits)
    if orbs == nothing
      G = Omega.group
      orbs = []
      for p in Omega.seeds
        if all(o -> !(p in o), orbs)
          push!(orbs, Hecke.orbit(Omega, p))
        end
      end
      AbstractAlgebra.set_special(Omega, :orbits => orbs)
    end
    return orbs
end


#############################################################################
##
##  `:elements` an array of points

function elements(Omega::GSetByElements)
    elms = AbstractAlgebra.get_special(Omega, :elements)
    if elms == nothing
      orbs = orbits(Omega)
      elms = union(map(elements, orbs)...)
      AbstractAlgebra.set_special(Omega, :elements => elms)
    end
    return elms
end


#############################################################################
##

# In fact, '<:GAPGroup' is not used at all in this function.
function permutation(Omega::GSetByElements{T}, g::BasicGAPGroupElem{T}) where T<:GAPGroup
    omega_list = GAP.julia_to_gap(elements(Omega))
    gfun = GAP.julia_to_gap(Omega.action_function)

    # The following works only because GAP does not check
    # whether the given group element 'g' is a group element.
    pi = GAP.Globals.PermutationOp(g, omega_list, gfun)

    sym = AbstractAlgebra.get_special(Omega, :action_range)
    if sym == nothing
      sym = symmetric_group(length(Omega))
      AbstractAlgebra.set_special(Omega, :action_range => sym )
    end
    return group_element(sym, pi)
end


############################################################################
##
##  action homomorphisms

# The following must be executed at runtime.
# (Eventually put it into some `__init__`.)
function __init_JuliaData()
    if ! GAP.Globals.ISB_GVAR(GAP.julia_to_gap("JuliaData"))
      GAP.evalstr("""
DeclareAttribute( "JuliaData", IsObject );");
InstallOtherMethod( ImagesRepresentative,
[ IsActionHomomorphism and HasJuliaData, IsMultiplicativeElementWithInverse ],
function( hom, elm )
local data;
data:= JuliaData( hom );
return Julia.Oscar.permutation(data[1], Julia.Oscar.group_element(data[2], elm)).X;
end );
""")
    end
end

function action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup
    acthom = AbstractAlgebra.get_special(Omega, :action_homomorphism)
    if acthom == nothing
      G = Omega.group
      omega_list = GAP.julia_to_gap(elements(Omega))
      gap_gens = map(x -> x.X, gens(G))
      gfun = GAP.julia_to_gap(Omega.action_function)

      # The following works only because GAP does not check
      # whether the given generators in GAP and Julia fit together.
      acthom = GAP.Globals.ActionHomomorphism(G.X, omega_list, GAP.julia_to_gap(gap_gens), GAP.julia_to_gap(gens(G)), gfun)

      # The first difficulty on the GAP side is `ImagesRepresentative`
      # (which is the easy direction of the action homomorphism):
      # In the method in question, GAP does not really know how to compute
      # the group element that actually acts from the given group element;
      # there is only a rudimentary `FunctionAction` inside the
      # `UnderlyingExternalSet` of the GAP homomorphism object `acthom`.
      # We could replace this function here,
      # but this would introduce overhead for mapping each point.
      # Thus we install a special `ImagesRepresentative` method in GAP;
      # note that we know how to get the Julia "actor" from the GAP group
      # element, by wrapping it into the corresponding Julia group element.
      # (Yes, this is also overhead.
      # The alternative would be to create a new type of Oscar homomorphism,
      # which uses `permutation` or something better for mapping elements.)
      __init_JuliaData()
      GAP.Globals.SetJuliaData(acthom, GAP.julia_to_gap([Omega, G]))

      sym = AbstractAlgebra.get_special(Omega, :action_range)
      if sym == nothing
        sym = symmetric_group(length(Omega))
        AbstractAlgebra.set_special(Omega, :action_range => sym )
      end
      acthom = GAPGroupHomomorphism{T,PermGroup}(G, sym, acthom)
      AbstractAlgebra.set_special(Omega, :action_homomorphism => acthom )
    end
    return acthom
end


"""
    isconjugate(Omega::GSet, omega1, omega2)

Return `true` if `omega1`, `omega2` are in the same orbit of `Omega`,
and `false` otherwise.
"""
isconjugate(Omega::GSet, omega1, omega2) = omega2 in orbit(Omega, omega1)


"""
    representative_action(Omega::GSet, omega1, omega2)

Determine whether `omega1`, `omega2` are in the same orbit of `Omega`.
If yes, return `true, g` where `g` is an element in the group `G` of
`Omega` that maps `omega1` to `omega2`.
If not, return `false, nothing`.
"""
function representative_action(Omega::GSet, omega1, omega2)
    # We do not call GAP's 'RepresentativeAction' with points, generators,
    # and actors.
    # The method in question would create a new 'ExternalSet' object
    # with a useless 'FunctionAction' value.
    # Instead, we delegate to the image of the action homomorphism.
    # (For that, we write down the elements of the G-set.
    # Computing the orbit of `omega1` or `omega2` would in principle suffice.)
    G = Omega.group
    acthom = action_homomorphism(Omega)
    elms = elements(Omega)
    pos1 = findfirst(isequal(omega1), elms)
    pos1 == nothing && return false, one(G)
    pos2 = findfirst(isequal(omega2), elms)
    pos2 == nothing && return false, one(G)
    img = GAP.Globals.RepresentativeAction(image(acthom)[1].X, pos1, pos2)
    img == GAP.Globals.fail && return false, one(G)
    pre = haspreimage(acthom, group_element(image(acthom)[1], img))
    @assert(pre[1])
    return true, pre[2]
end


############################################################################

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

