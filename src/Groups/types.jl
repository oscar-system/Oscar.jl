import Hecke:
    abelian_group,
    automorphism_group,
    center,
    codomain,
    cokernel,
    compose,
    degree,
    derived_series,
    direct_product,
    domain,
    elem_type,
    elements,
    free_abelian_group,
    gen,
    gens,
    haspreimage,
    hom,
    id_hom,
    image,
    index,
    inv!,
    isabelian,
    isbijective,
    ischaracteristic,
    isconjugate,
    iscyclic,
    isinjective,
    isinvertible,
    isisomorphic,
    isnormal,
    issimple,
    issubgroup,
    issurjective,
    kernel,
    Map,
    mul,
    mul!,
    ngens,
    normal_closure,
    one!,
    order,
    parent_type,
    perm,
    preimage,
    quo,
    representative,
    SetMap,
    small_group,
    sub,
    subgroups

import Base: ==, parent, show

import GAP.GapObj

export
    AutomorphismGroup,
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

"""
TODO: document this
"""
abstract type GAPGroup <: AbstractAlgebra.Group end
abstract type GAPGroupElem{T<:GAPGroup} <: AbstractAlgebra.GroupElem end

"""
    BasicGAPGroupElem{T<:GAPGroup}
The type `BasicGAPGroupElem` gathers all types of group elements described only by an underlying GAP object.
"""
struct BasicGAPGroupElem{T<:GAPGroup} <: GAPGroupElem{T}
   parent::T
   X::GapObj
end

Base.hash(x::GAPGroup, h::UInt) = h # FIXME
Base.hash(x::GAPGroupElem, h::UInt) = h # FIXME


"""
    PermGroup

Groups of permutations. Every group of this type is the subgroup of Sym(n) for some n.

# Examples
- `symmetric_group(n::Int)`: the symmetric group Sym(n)
- `alternating_group(n::Int)`: the alternating group Alt(n)
- subgroups of Sym(n)
- `dihedral_group(PermGroup, n::Int)`: the dihedral group D(n) as group of permutations. Same holds replacing `dihedral_group` by `quaternion_group`
"""
mutable struct PermGroup <: GAPGroup
   X::GapObj
   deg::Int64       # G < Sym(deg)
   AbstractAlgebra.@declare_other
   
   function PermGroup(G::GapObj)
     @assert GAP.Globals.IsPermGroup(G)
     n = GAP.gap_to_julia(Int64, GAP.Globals.LargestMovedPoint(G))
     z = new(G, n)
     return z
   end
   
   function PermGroup(G::GapObj, deg::Int)
     @assert GAP.Globals.IsPermGroup(G)
     z = new(G, deg)
     return z
   end
end

"""
    PermGroupElem

Element of a group of permutation. It is displayed as product of disjoint cycles.
# Assumptions:
- for `x`,`y` in Sym(n), the product `xy` is read from left to right;
- for `x` in Sym(n) and `i` in {1,...,n}, `i^x` and `x(i)` return the image of `i` under the action of `x`.
"""
const PermGroupElem = BasicGAPGroupElem{PermGroup}

"""
    PcGroup

Polycyclic group
# Examples:
- `cyclic_group(n::Int)`: cyclic group of order `n`
- `abelian_group(v::Vector{Int})`: direct product of cyclic groups of order v[1],v[2],...,v[length(v)]
"""
mutable struct PcGroup <: GAPGroup
  X::GapObj
  AbstractAlgebra.@declare_other

  function PcGroup(G::GapObj)
    @assert GAP.Globals.IsPcGroup(G)
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

Finitely presented group. It can be defined via the function ``free_group``.
"""
mutable struct FPGroup <: GAPGroup
  X::GapObj
  AbstractAlgebra.@declare_other
  
  function FPGroup(G::GapObj)
    @assert GAP.Globals.IsFpGroup(G)
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
#  Group Homomorphism
#
################################################################################

abstract type GAPMap <: SetMap end

struct GAPGroupHomomorphism{S<: GAPGroup, T<: GAPGroup} <: Map{S,T,GAPMap,GAPGroupHomomorphism{S,T}}
   domain::S
   codomain::T
   map::GapObj
end

"""
    AutomorphismGroup{T} <: GAPGroup

Group of automorphisms over a group of type `T`. It can be defined via the function ``automorphism_group``.
"""
mutable struct AutomorphismGroup{T} <: GAPGroup
  X::GapObj
  G::T
  AbstractAlgebra.@declare_other

  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAP.Globals.IsGroupOfAutomorphisms(G)
    z = new{T}(G, H)
    return z
  end
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

Semidirect product of two groups of type `S` and `T` respectively, or subgroup of a semidirect product of groups.
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

Wreath product of a group `G` and a group of permutations `H`, or a generic group `H` together with the homomorphism `a` from `H` to a permutation group.
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
TODO: document this

`elem_type` maps a group to the type of its elements. For now,
a group of type `T` has elements of type `GAPGroupElem{T}`. So
we provide it mostly for consistency with other parts of OSCAR.
In the future, a more elaborate setup for group element types
might also be needed.
"""

elem_type(::Type{T}) where T <: GAPGroup = BasicGAPGroupElem{T}



#
# The array _gap_group_types contains pairs (X,Y) where
# X is a GAP filter such as IsPermGroup, and Y is a corresponding
# Julia type such as `PermGroup`.
#
 
const _gap_group_types = []

function _get_type(G::GapObj)
  for pair in _gap_group_types
    if pair[1](G)
      return pair[2]
    end
  end
  error("Not a known type of group")
end


