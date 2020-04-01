export PermGroup, PermGroupElem,
       MatrixGroup, MatrixGroupElem,
       PcGroup, PcGroupElem,
       FPGroup, FPGroupElem,
       AutomorphismGroup, AutomorphismGroupElem,
       elem_type

"""
TODO: document this
"""
abstract type Group <: AbstractAlgebra.Group end

"""
TODO: document this
"""
struct GroupElem{T<:Group} <: AbstractAlgebra.GroupElem
   parent::T
   X::GapObj
end

"""
TODO: document this
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
const PermGroupElem = GroupElem{PermGroup}

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
const MatrixGroupElem = GroupElem{MatrixGroup}

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
const PcGroupElem = GroupElem{PcGroup}

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
const FPGroupElem = GroupElem{FPGroup}

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
const AutomorphismGroupElem = GroupElem{AutomorphismGroup}


"""
TODO: document this

`elem_type` maps a group to the type of its elements. For now,
a group of type `T` has elements of type `GroupElem{T}`. So
we provide it mostly for consistency with other parts of OSCAR.
In the future, a more elaborate setup for group element types
might also be needed.
"""

elem_type(::T) where T <: Group = GroupElem{T}


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
