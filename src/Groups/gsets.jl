# This file contains code related to G-sets
# The idea of the implementation is that available GAP functionality
# for computing orbits, permutations, actions, etc. can be used,
# but that local replacements by pure Julia code (independent of GAP)
# are welcome.

import Hecke.orbit

export
    all_blocks,
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

@attributes mutable struct GSetByElements{T} <: GSet{T}
    group::T
    action_function::Function
    seeds

    function GSetByElements(G::T, fun::Function, seeds) where T<:GAPGroup
        @assert ! isempty(seeds)
        return new{T}(G, fun, seeds, Dict{Symbol,Any}())
    end
end
#TODO: How can I specify that `seeds` should be an iterable object?

# TODO: document `acting_group` and GSets in general
acting_group(gset::GSetByElements) = gset.group


#############################################################################
##
##  general method with explicit action function

"""
    gset(G::GAPGroup[, fun::Function], Omega)

Return the G-set that consists of the closure of the seeds `Omega`
under the action of `G` defined by `fun`.

This means that the result contains all elements `fun(omega, g)`
for `omega` in `Omega` and `g` in `G`.

`fun` can be omitted if the element type of `Omega` implies
a reasonable default,
for example, if `G` is a `PermGroup` and `Omega` is a `Vector{T}` where
`T` is one of `Int`, `Set{Int}`, `Vector{Int}`.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> length(gset(G, [[1]]))  # natural action
4

julia> length(gset(G, [[1, 2]]))  # action on ordered pairs
12

julia> length(gset(G, on_sets, [[1, 2]]))  # action on unordered pairs
6

```
"""
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
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:IntegerUnion = GSetByElements(G, ^, Omega)

## action of permutations on sets of positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Set{T2} where T2<:IntegerUnion = GSetByElements(G, on_sets, Omega)

## action of permutations on vectors of positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Vector{T2} where T2<:IntegerUnion = GSetByElements(G, on_tuples, Omega)

## action of permutations on tuples of positive integers
gset_by_type(G::PermGroup, Omega, ::Type{T}) where T<:Tuple{Vararg{T2}} where T2<:IntegerUnion = GSetByElements(G, on_tuples, Omega)

## (add more such actions: on sets of sets, on sets of tuples, ...)

## natural action of a permutation group on the integers 1, ..., degree
function gset(G::PermGroup)
    omega = gset(G, 1:G.deg)
    set_attribute!(omega, :elements => omega.seeds)
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
    set_attribute!(omega, :elements => omega.seeds)
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

orbit(omega::ElementOfGSet) = orbit(omega.gset, omega.obj)


unwrap(omega::Any) = omega

unwrap(omega::ElementOfGSet) = omega.obj


#############################################################################
##
##  `:orbit`

"""
    orbit(G::PermGroup[, fun::Function], omega)

Return the G-set that consists of the images of `omega`
under the action of `G` defined by `fun`.

This means that the result contains all elements `fun(omega, g)`
for `g` in `G`.

`fun` can be omitted if the type of `Omega` implies a reasonable default,
for example, if `G` is a `PermGroup` and `omega` is
one of `Int`, `Set{Int}`, `Vector{Int}`.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> length(orbit(G, 1))
4

julia> length(orbit(G, [1, 2]))
12

julia> length(orbit(G, on_sets, [1, 2]))
6

```
"""
orbit(G::PermGroup, omega) = gset_by_type(G, [omega], typeof(omega))

orbit(G::PermGroup, fun::Function, omega) = GSetByElements(G, fun, [omega])

"""
    orbit(Omega::GSet, omega::T) where T

Return the G-set that consists of the elements `fun(omega, g)` where
`g` is in the group of `Omega` and `fun` is the underlying action of `Omega`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4), (5,6) ])

julia> Omega = gset(G, [1, 5]);

julia> length(orbit(Omega, 1))
4

```
"""
function orbit(Omega::GSetByElements{<:GAPGroup}, omega::T) where T
    G = Omega.group
    acts = GapObj(gens(G))
    gfun = GapObj(Omega.action_function)

    # The following works only because GAP does not check
    # whether the given (dummy) group 'G.X' fits to the given generators,
    # or whether the elements of 'acts' are group elements.
    orb = Vector{T}(GAP.Globals.Orbit(G.X, omega, acts, acts, gfun)::GapObj)

    res = as_gset(Omega.group, Omega.action_function, orb)
    # We know that this G-set is transitive.
    set_attribute!(res, :orbits => [orb])
    return res
end
#T check whether omega lies in Omega?

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
    set_attribute!(res, :orbits => [orbarray])
    return res
end


#############################################################################
##
##  `:orbits` a vector of G-sets

"""
    orbits(Omega::GSet)

Return the vector of transitive G-sets in `Omega`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4), (5,6) ])

julia> orbs = orbits(gset(G));

julia> map(elements, orbs)
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [5, 6]

```
"""
@attr Vector{GSetByElements{TG}} function orbits(Omega::T) where T <: GSetByElements{TG} where TG <: GAPGroup
  G = Omega.group
  orbs = T[]
  for p in Omega.seeds
    if all(o -> !(p in o), orbs)
      push!(orbs, orbit(Omega, p))
    end
  end
  return orbs
end

"""
    orbits(G::PermGroup)

Return the orbits of the natural G-set of `G`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4), (5,6) ])

julia> orbs = orbits(G);

julia> map(length, orbs)
2-element Vector{Int64}:
 4
 2

```
"""
@attr Vector{GSetByElements{PermGroup}} orbits(G::PermGroup) = orbits(gset(G))


#############################################################################
##
##  `:elements` a vector of points

@attr Any function elements(Omega::GSetByElements)
  orbs = orbits(Omega)
  elms = union(map(elements, orbs)...)
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

    sym = get_attribute!(Omega, :action_range) do
      return symmetric_group(length(Omega))
    end
    return group_element(sym, pi)
end


############################################################################
##
##  action homomorphisms

# Use a GAP attribute for caching the mapping.
# The following must be executed at runtime,
# the function gets called in Oscar's `__init__`.
function __init_JuliaData()
    if ! hasproperty(GAP.Globals, :JuliaData)
      GAP.evalstr("""
DeclareAttribute( "JuliaData", IsObject );
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

"""
    action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup

Return the group homomorphism `act` with domain `G = Omega.group`
and codomain `symmetric_group(n)` that describes the action of `G` on `Omega`,
where `Omega` has `n` elements.

This means that if an element `g` in `G` maps `elements(Omega)[i]` to
`elements(Omega)[j]` then `act(g)` maps `i` to `j`.

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> Omega = gset(G, [Set([1, 2])]);  # action on unordered pairs

julia> acthom = action_homomorphism(Omega)
Group homomorphism from 
Sym( [ 1 .. 6 ] )
to
Sym( [ 1 .. 15 ] )

julia> g = gens(G)[1]
(1,2,3,4,5,6)

julia> elms = elements(Omega);

julia> actg = acthom(g)
(1,2,3,5,7,10)(4,6,8,11,14,13)(9,12,15)

julia> elms[1]^g == elms[2]
true

julia> 1^actg == 2
true

```
"""
@attr GAPGroupHomomorphism{T, PermGroup} function action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup
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
  GAP.Globals.SetJuliaData(acthom, GAP.julia_to_gap([Omega, G]))

  sym = get_attribute!(Omega, :action_range) do
    return symmetric_group(length(Omega))
  end
  return GAPGroupHomomorphism(G, sym, acthom)
end

# for convenience: create the G-set on the fly
function action_homomorphism(G::PermGroup, Omega)
  return action_homomorphism(gset_by_type(G, Omega, eltype(Omega)))
end

function action_homomorphism(G::PermGroup, fun::Function, Omega)
  return action_homomorphism(GSetByElements(G, fun, Omega))
end


"""
    isconjugate(Omega::GSet, omega1, omega2)

Return `true` if `omega1`, `omega2` are in the same orbit of `Omega`,
and `false` otherwise.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4), (5,6) ])

julia> Omega = gset(G);

julia> isconjugate(Omega, 1, 2)
true

julia> isconjugate(Omega, 1, 5)
false

```
"""
isconjugate(Omega::GSet, omega1, omega2) = omega2 in orbit(Omega, omega1)


"""
    representative_action(Omega::GSet, omega1, omega2)

Determine whether `omega1`, `omega2` are in the same orbit of `Omega`.
If yes, return `true, g` where `g` is an element in the group `G` of
`Omega` that maps `omega1` to `omega2`.
If not, return `false, nothing`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4), (5,6) ])

julia> Omega = gset(G);

julia> representative_action(Omega, 1, 2)
(true, (1,2))

julia> representative_action(Omega, 1, 5)
(false, ())

```
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

Base.getindex(Omega::GSet, i::Int) = elements(Omega)[i]

# TODO:
# - Base.iterate


blocks(G::GSet) = error("not implemented")
maximal_blocks(G::GSet) = error("not implemented")
representatives_minimal_blocks(G::GSet) = error("not implemented")
all_blocks(G::GSet) = error("not implemented")

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
    blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))

Return a G-set that is a block system for the action of `G` on `L`,
i.e., a non-trivial partition of `L` preserved by the action of `G`.

Here, `L` must be a subvector of `1:degree(G)` on which `G` acts transitively.
`G` may move points outside `L`, in this case the restriction of the action
of the set stabilizer of `L` in `G` to `L` is considered.

An exception is thrown if this action is not transitive.

# Examples
```jldoctest
julia> g = sylow_subgroup(symmetric_group(4), 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4) ])

julia> elements(blocks(g))
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3, 4]

```
"""
function blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
   @assert istransitive(G, L) "The group action is not transitive"
   bl = Vector{Vector{Int}}(GAP.Globals.Blocks(G.X, GapObj(L))::GapObj)
   omega = gset(G, bl)
   set_attribute!(omega, :elements => omega.seeds)
   return omega
end

"""
    maximal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))

Return a G-set that is a maximal block system for the action of `G` on `L`,
i.e., a maximal non-trivial partition of `L` preserved by the action of `G`.

Here, `L` must be a subvector of `1:degree(G)` on which `G` acts transitively.
`G` may move points outside `L`, in this case the restriction of the action
of the set stabilizer of `L` in `G` to `L` is considered.

An exception is thrown if this action is not transitive.

# Examples
```jldoctest
julia> G = transitive_group(8, 2)
4[x]2

julia> elements(maximal_blocks(G))
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 8]
 [4, 5, 6, 7]

```
"""
function maximal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
   @assert istransitive(G, L) "The group action is not transitive"
   bl = Vector{Vector{Int}}(GAP.Globals.MaximalBlocks(G.X, GapObj(L))::GapObj)
   omega = gset(G, bl)
   set_attribute!(omega, :elements => omega.seeds)
   return omega
end


"""
    representatives_minimal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))

Return a vector of block representatives for all minimal non-trivial block
systems for the action of `G` on `L`.

Here, `L` must be a subvector of `1:degree(G)` on which `G` acts transitively.
`G` may move points outside `L`, in this case the restriction of the action
of the set stabilizer of `L` in `G` to `L` is considered.

An exception is thrown if this action is not transitive.

# Examples
```jldoctest
julia> G = transitive_group(8, 2)
4[x]2

julia> representatives_minimal_blocks(G)
3-element Vector{Vector{Int64}}:
 [1, 3]
 [1, 5]
 [1, 7]

```
"""
function representatives_minimal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
   @assert istransitive(G, L) "The group action is not transitive"
   return Vector{Vector{Int}}(GAP.Globals.RepresentativesMinimalBlocks(G.X, GapObj(L))::GapObj)
end


"""
    all_blocks(G::PermGroup)

Return a vector of representatives of all block systems for the action of
`G` on the set of moved points of `G`.

# Examples
```jldoctest
julia> G = transitive_group(8, 2)
4[x]2

julia> all_blocks(G)
6-element Vector{Vector{Int64}}:
 [1, 2, 3, 8]
 [1, 5]
 [1, 3, 5, 7]
 [1, 3]
 [1, 3, 4, 6]
 [1, 7]

```
"""
all_blocks(G::PermGroup) = Vector{Vector{Int}}(GAP.Globals.AllBlocks(G.X))
#TODO: Do we really want to act on the set of moved points?


"""
    transitivity(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return the maximum `k` such that the action of `G` on `L` is
`k`-transitive. The output is `0` if `G` is not transitive on `L`.

# Examples
```jldoctest
julia> transitivity(mathieu_group(24))
5

julia> transitivity(symmetric_group(6))
6

```
"""
transitivity(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAP.Globals.Transitivity(G.X, GapObj(L))::Int

"""
    istransitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is transitive.

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> istransitive(G)
true

julia> istransitive(sylow_subgroup(G, 2)[1])
false

```
"""
istransitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsTransitive(G.X, GapObj(L))
# Note that this definition does not coincide with that of the
# property `GAP.Globals.IsTransitive`, for which the default domain
# of the action is the set of moved points.


"""
    isprimitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is primitive, that is,
the action is transitive and the point stabilizers are maximal in `G`.

# Examples
```jldoctest
julia> G = alternating_group(6);

julia> mx = filter(istransitive, maximal_subgroup_reps(G))
3-element Vector{PermGroup}:
 Group([ (1,2)(3,4), (1,2)(5,6), (1,3,5)(2,4,6), (1,3)(2,4) ])
 Group([ (1,2,3), (4,5,6), (1,2)(4,5), (1,5,2,4)(3,6) ])
 PSL(2,5)

julia> [(order(H), isprimitive(H)) for H in mx]
3-element Vector{Tuple{fmpz, Bool}}:
 (24, 0)
 (36, 0)
 (60, 1)

```
"""
isprimitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsPrimitive(G.X, GapObj(L))


"""
    isregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is regular
(i.e., transitive and semiregular).

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> H = sub(G, [G([2, 3, 4, 5, 6, 1])])[1]
Group([ (1,2,3,4,5,6) ])

julia> isregular(H)
true

julia> isregular(G)
false

```
"""
isregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsRegular(G.X, GapObj(L))


"""
    issemiregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is semiregular
(i.e., the stabilizer of each point is the identity).

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> H = sub(G, [G([2, 3, 1, 5, 6, 4])])[1]
Group([ (1,2,3)(4,5,6) ])

julia> issemiregular(H)
true

julia> isregular(H)
false

```
"""
issemiregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsSemiRegular(G.X, GapObj(L))
