export
    AutomorphismGroup,
    AutomorphismGroupElem,
    DirectProductGroup,
    DirectProductOfElem,
    FPGroup,
    FPGroupElem,
    GAPGroupElem,
    GAPGroupHomomorphism,
    PcGroup,
    PcGroupElem,
    PermGroup,
    PermGroupElem,
    SemidirectProductGroup,
    WreathProductGroup

@doc Markdown.doc"""
    GAPGroup <: AbstractAlgebra.Group

Each object of the abstract type `GAPGroup` stores a group object from
the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `G` is a `GAPGroup`, then `GapObj(G)` is the `GapObj` underlying `G`.

Concrete subtypes of `GAPGroup` are `PermGroup`, `FPGroup`, `PcGroup`,
and `MatrixGroup`.
"""
abstract type GAPGroup <: AbstractAlgebra.Group end

## `GapGroup` to GAP group
GAP.julia_to_gap(obj::GAPGroup) = obj.X

@doc Markdown.doc"""
    GAPGroupElem <: AbstractAlgebra.GroupElem

Each object of the abstract type `GAPGroupElem` stores a group element
object from the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `g` is a `GAPGroupElem`, then `GapObj(g)` is the `GapObj` underlying `g`.
"""
abstract type GAPGroupElem{T<:GAPGroup} <: AbstractAlgebra.GroupElem end

## `GapGroupElem` to GAP group element
GAP.julia_to_gap(obj::GAPGroupElem) = obj.X

@doc Markdown.doc"""
    BasicGAPGroupElem{T<:GAPGroup} <: GAPGroupElem{T}

The type `BasicGAPGroupElem` gathers all types of group elements
described *only* by an underlying GAP object.

If $x$ is an element of the group `G` of type `T`,
then the type of $x$ is `BasicGAPGroupElem{T}`.
"""
struct BasicGAPGroupElem{T<:GAPGroup} <: GAPGroupElem{T}
   parent::T
   X::GapObj
end

function Base.deepcopy_internal(x::BasicGAPGroupElem, dict::IdDict)
  X = Base.deepcopy_internal(x.X, dict)
  return BasicGAPGroupElem(x.parent, X)
end

Base.hash(x::GAPGroup, h::UInt) = h # FIXME
Base.hash(x::GAPGroupElem, h::UInt) = h # FIXME


"""
    PermGroup

Groups of permutations.
Every group of this type is the subgroup of Sym(n) for some n.

# Examples
- `symmetric_group(n::Int)`: the symmetric group Sym(n)
- `alternating_group(n::Int)`: the alternating group Alt(n)
- subgroups of Sym(n)
- `dihedral_group(PermGroup, n::Int)`:
  the dihedral group of order `n` as a group of permutations.
  Same holds replacing `dihedral_group` by `quaternion_group`
"""
@attributes mutable struct PermGroup <: GAPGroup
   X::GapObj
   deg::Int64       # G < Sym(deg)
   
   function PermGroup(G::GapObj)
     @assert GAPWrap.IsPermGroup(G)
     n = GAPWrap.LargestMovedPoint(G)::Int
     if n == 0
       # We support only positive degrees.
       # (`symmetric_group(0)` yields an error,
       # and `symmetric_group(1)` yields a GAP group with `n == 0`.)
       n = 1
     end
     z = new(G, n)
     return z
   end
   
   function PermGroup(G::GapObj, deg::Int)
     @assert GAPWrap.IsPermGroup(G) && deg > 0 && deg >= GAPWrap.LargestMovedPoint(G)::Int
     z = new(G, deg)
     return z
   end
end

"""
    PermGroupElem

Element of a group of permutations.
It is displayed as product of disjoint cycles.
# Assumptions:
- for `x`,`y` in Sym(n), the product `xy` is read from left to right;
- for `x` in Sym(n) and `i` in {1,...,n}, `i^x` and `x(i)` return the image of `i` under the action of `x`.
"""
const PermGroupElem = BasicGAPGroupElem{PermGroup}

"""
    PcGroup

Polycyclic group

# Examples
- `cyclic_group(n::Int)`: cyclic group of order `n`
- `abelian_group(PcGroup, v::Vector{Int})`:
  direct product of cyclic groups of the orders
  `v[1]`, `v[2]`, ..., `v[length(v)]`
"""
@attributes mutable struct PcGroup <: GAPGroup
  X::GapObj

  function PcGroup(G::GapObj)
    @assert GAPWrap.IsPcGroup(G)
    z = new(G)
    return z
  end
end

"""
    PcGroupElem

Element of a polycyclic group.
"""
const PcGroupElem = BasicGAPGroupElem{PcGroup}

"""
    FPGroup

Finitely presented group.
Such groups can be constructed a factors of free groups,
see [`free_group`](@ref).
"""
@attributes mutable struct FPGroup <: GAPGroup
  X::GapObj
  
  function FPGroup(G::GapObj)
    @assert GAPWrap.IsSubgroupFpGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
const FPGroupElem = BasicGAPGroupElem{FPGroup}


################################################################################
#
# Construct an Oscar group wrapping the GAP group `obj`
# *and* compatible with a given Oscar group `G`.

# default: ignore `G`
_oscar_group(obj::GapObj, G::T) where T <: GAPGroup = T(obj)

# `PermGroup`: set the degree of `G`
function _oscar_group(obj::GapObj, G::PermGroup)
  n = GAP.Globals.LargestMovedPoint(obj)
  N = degree(G)
  n <= N || error("requested degree ($N) is smaller than the largest moved point ($n)")
  return PermGroup(obj, N)
end

# `MatrixGroup`: set dimension and ring of `G`
# (This cannot be defined here because `MatrixGroup` is not yet defined.)


################################################################################
#
# "Coerce" an Oscar group `G` to one that is compatible with
# the given Oscar group `S`.
compatible_group(G::T, S::T) where T <: GAPGroup = _oscar_group(G.X, S)


################################################################################
#
#  Group Homomorphism
#
################################################################################

abstract type GAPMap <: SetMap end

struct GAPGroupHomomorphism{S<: GAPGroup, T<: GAPGroup} <: Map{S,T,GAPMap,GAPGroupHomomorphism{S,T}}
   domain::S
   codomain::T
   map::GapObj

   function GAPGroupHomomorphism(G::S, H::T, mp::GapObj) where {S<: GAPGroup, T<: GAPGroup}
     return new{S, T}(G, H, mp)
   end
end


"""
    AutomorphismGroup{T} <: GAPGroup

Group of automorphisms over a group of type `T`. It can be defined via the function `automorphism_group`.
"""
@attributes mutable struct AutomorphismGroup{T} <: GAPGroup
  X::GapObj
  G::T

  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAPWrap.IsGroupOfAutomorphisms(G)
    z = new{T}(G, H)
    return z
  end
end

function AutomorphismGroup(G::GapObj, H::T) where T
  return AutomorphismGroup{T}(G, H)
end

(aut::AutomorphismGroup{T} where T)(x::GapObj) = group_element(aut,x)

const AutomorphismGroupElem{T} = BasicGAPGroupElem{AutomorphismGroup{T}} where T

function Base.show(io::IO, AGE::AutomorphismGroupElem{GrpAbFinGen}) 
    println(io, "Automorphism of ", GrpAbFinGen, " with matrix representation ", matrix(AGE))
end


################################################################################
#
#  Composite Groups
#
################################################################################


"""
    DirectProductGroup

Either direct product of two or more groups of any type, or subgroup of a direct product of groups.
"""
struct DirectProductGroup <: GAPGroup
  X::GapObj
  L::Vector{<:GAPGroup}   # list of groups
  Xfull::GapObj      # direct product of the GAP groups of L
  isfull::Bool     # true if G is direct product of the groups of L, false if it is a proper subgroup
end


"""
    SemidirectProductGroup{S,T}

Semidirect product of two groups of type `S` and `T` respectively, or
subgroup of a semidirect product of groups.
"""
struct SemidirectProductGroup{S<:GAPGroup, T<:GAPGroup} <: GAPGroup 
  X::GapObj
  N::S              # normal subgroup
  H::T              # group acting on N
  f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}        # action of H on N
  Xfull::GapObj         # full semidirect product: X is a subgroup of Xfull. 
  isfull::Bool     # true if X==Xfull
end

"""
    WreathProductGroup

Wreath product of a group `G` and a group of permutations `H`, or a generic
group `H` together with the homomorphism `a` from `H` to a permutation
group.
"""
struct WreathProductGroup <: GAPGroup
  X::GapObj
  G::GAPGroup
  H::GAPGroup
  a::GAPGroupHomomorphism   # morphism from H to the permutation group
  Xfull::GapObj            # if H does not move all the points, this is the wreath product of (G, Sym(degree(H))
  isfull::Bool             # true if Xfull == X
end


"""
    elem_type(::Type{T}) where T <: GAPGroup
    elem_type(::T) where T <: GAPGroup

`elem_type` maps (the type of) a group to the type of its elements.
For now, a group of type `T` has elements of type `BasicGAPGroupElem{T}`.
So we provide it mostly for consistency with other parts of OSCAR.
In the future, a more elaborate setup for group element types
might also be needed.
"""
elem_type(::Type{T}) where T <: GAPGroup = BasicGAPGroupElem{T}
elem_type(::T) where T <: GAPGroup = BasicGAPGroupElem{T}

Base.eltype(::Type{T}) where T <: GAPGroup = BasicGAPGroupElem{T}

# `parent_type` is defined and documented in AbstractAlgebra.
parent_type(::Type{T}) where T<:BasicGAPGroupElem{S} where S = S
parent_type(::T) where T<:BasicGAPGroupElem{S} where S = S

#
# The array _gap_group_types contains pairs (X,Y) where
# X is a GAP filter such as IsPermGroup, and Y is a corresponding
# Julia type such as `PermGroup`.
#
 
const _gap_group_types = Tuple{GAP.GapObj, Type}[]

function _get_type(G::GapObj)
  for pair in _gap_group_types
    if pair[1](G)
      if pair[2] == MatrixGroup
#T HACK: We need more information in the case of matrix groups.
#T (Usually we should not need to guess the Oscar side of a GAP group.)
        return function(dom::GAP.GapObj)
                 deg = GAP.Globals.DimensionOfMatrixGroup(dom)
                 iso = iso_gap_oscar(GAP.Globals.FieldOfMatrixGroup(dom))
                 ring = codomain(iso)
                 matgrp = MatrixGroup(deg, ring)
                 matgrp.ring_iso = inv(iso)
                 matgrp.X = dom
                 return matgrp
               end
      else
        return pair[2]
      end
    end
  end
  error("Not a known type of group")
end
