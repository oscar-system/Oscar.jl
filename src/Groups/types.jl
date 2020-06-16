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
    mul,
    mul!,
    ngens,
    one!,
    order,
    perm,
    preimage,
    quo,
    representative,
    small_group,
    sub,
    subgroups

import Base: ==, parent, show

import GAP.GapObj

export
    AutomorphismGroup,
    DirectProductOfGroups,
    DirectProductOfElem,
    elem_type,
    FPGroup,
    FPGroupElem,
    MatrixGroup,
    MatrixGroupElem,
    PcGroup,
    PcGroupElem,
    PermGroup,
    PermGroupElem

"""
TODO: document this
"""
abstract type GAPGroup <: AbstractAlgebra.Group end
#abstract type GroupElem <: AbstractAlgebra.GroupElem

"""
TODO: document this
"""
struct GAPGroupElem{T<:GAPGroup} <: AbstractAlgebra.GroupElem
   parent::T
   X::GapObj
end

Base.hash(x::GAPGroup) = 0 # FIXME
Base.hash(x::GAPGroupElem) = 0 # FIXME


"""
    PermGroup
Groups of permutations. Every group of this type is the subgroup of Sym(n) for some n.

# Examples
- `symmetric_group(n::Int)`: the symmetric group Sym(n)
- `alternating_group(n::Int)`: the alternating group Alt(n)
- subgroups of Sym(n)
- `dihedral_group(PermGroup, n::Int)`: the dihedral group D(n) as group of permutations. Same holds replacing `dihedral_group` by `quaternion_group`
"""
struct PermGroup <: GAPGroup
   X::GapObj
   deg::Int64       # G < Sym(deg)
   
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
- for `x` in Sym(n) and `i` in {1,...,n}, `x(i)` return the image of `i` under the action of `x`.
"""
const PermGroupElem = GAPGroupElem{PermGroup}

"""
    MatrixGroup
Groups of matrices. Every group of this type is the subgroup of GL(n,q) for some integer `n` and prime power `q`

# Examples
- `GL(n::Int, q::Int)`: the general linear group GL(n,q)
- `SL(n::Int)`: the special linear group SL(n,q)
- groups of isometries
"""
struct MatrixGroup <: GAPGroup
  X::GapObj
  function MatrixGroup(G::GapObj)
    @assert GAP.Globals.IsMatrixGroup(G)
    z = new(G)
    return z
  end
end

"""
    MatrixGroupElem
Element of a matrix group.
"""
const MatrixGroupElem = GAPGroupElem{MatrixGroup}

#display(x::MatrixGroupElem) = GAP.Globals.Display(x.X)

"""
    PcGroup
Polycyclic group
# Examples:
- `cyclic_group(n::Int)`: cyclic group of order `n`
- `abelian_group(v::Vector{Int})`: direct product of cyclic groups of order v[1],v[2],...,v[length(v)]
"""
struct PcGroup <: GAPGroup
  X::GapObj
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
const PcGroupElem = GAPGroupElem{PcGroup}

"""
    FPGroup
Finitely presented group. It can be defined via the function ``free_group``.
"""
struct FPGroup <: GAPGroup
  X::GapObj
  
  function FPGroup(G::GapObj)
    @assert GAP.Globals.IsFpGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
const FPGroupElem = GAPGroupElem{FPGroup}

"""
    AutomorphismGroup{T} <: GAPGroup
Group of automorphisms over a group of type `T`. It can be defined via the function ``automorphism_group``.
"""
struct AutomorphismGroup{T} <: GAPGroup
  X::GapObj
  G::T
  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAP.Globals.IsGroupOfAutomorphisms(G)
    z = new{T}(G, H)
    return z
  end
end

"""
    DirectProductOfGroups{S,T}
Direct product of two groups of type `S` and `T` respectively, or subgroup of a direct product of groups.
"""
struct DirectProductOfGroups{S<:GAPGroup, T<:GAPGroup} <: GAPGroup
  X::GapObj
  G1::S
  G2::T
  IsFull::Bool     # true if it is G1xG2; false if it is a proper subgroup of G1xG2

  function DirectProductOfGroups(G::GapObj, G1::S, G2::T, isf::Bool) where S where T
    z = new{S,T}(G,G1,G2,isf)
    return z
  end
end


"""
TODO: document this

`elem_type` maps a group to the type of its elements. For now,
a group of type `T` has elements of type `GAPGroupElem{T}`. So
we provide it mostly for consistency with other parts of OSCAR.
In the future, a more elaborate setup for group element types
might also be needed.
"""

elem_type(::T) where T <: GAPGroup = GAPGroupElem{T}


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

################################################################################
#
#  Group Homomorphism
#
################################################################################

struct GAPGroupHomomorphism{S<: GAPGroup, T<: GAPGroup}
   domain::S
   codomain::T
   map::GapObj
end
