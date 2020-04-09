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


export PermGroup, PermGroupElem
export MatrixGroup, MatrixGroupElem
export PcGroup, PcGroupElem
export FPGroup, FPGroupElem
export AutomorphismGroup, AutomorphismGroupElem
export elem_type

"""
TODO: document this
"""
abstract type Group <: AbstractAlgebra.Group end

"""
TODO: document this
"""
struct GAPGroupElem{T<:Group} <: AbstractAlgebra.GroupElem
   parent::T
   X::GapObj
end

"""
    PermGroup
Groups of permutations. Every group of this type is the subgroup of Sym(n) for some n.
#Examples
- `symmetric_group(n::Int)`: the symmetric group Sym(n)
- `alternating_group(n::Int)`: the alternating group Alt(n)
- subgroups of Sym(n)
- `dihedral_group(PermGroup, n::Int)`: the dihedral group D(n) as group of permutations. Same holds replacing `dihedral_group` by `quaternion_group`
"""
struct PermGroup <: Group
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
TODO: document this
"""
const PermGroupElem = GAPGroupElem{PermGroup}

"""
TODO: document this
"""
struct MatrixGroup <: Group
  X::GapObj
  function MatrixGroup(G::GapObj)
    @assert GAP.Globals.IsMatrixGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
const MatrixGroupElem = GAPGroupElem{MatrixGroup}

"""
TODO: document this
"""
struct PcGroup <: Group
  X::GapObj
  function PcGroup(G::GapObj)
    @assert GAP.Globals.IsPcGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
const PcGroupElem = GAPGroupElem{PcGroup}

"""
TODO: document this
"""
struct FPGroup <: Group
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
TODO: document this
"""
struct AutomorphismGroup{T} <: Group
  X::GapObj
  G::T
  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAP.Globals.IsAutomorphismGroup(G)
    z = new{T}(G, H)
    return z
  end
end

"""
TODO: document this
"""
const AutomorphismGroupElem = GAPGroupElem{AutomorphismGroup}


"""
TODO: document this

`elem_type` maps a group to the type of its elements. For now,
a group of type `T` has elements of type `GAPGroupElem{T}`. So
we provide it mostly for consistency with other parts of OSCAR.
In the future, a more elaborate setup for group element types
might also be needed.
"""

elem_type(::T) where T <: Group = GAPGroupElem{T}


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

struct GAPGroupHomomorphism{S<: Group, T<: Group}
   domain::S
   codomain::T
   map::GapObj
end
