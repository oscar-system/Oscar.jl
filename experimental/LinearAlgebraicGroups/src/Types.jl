###############################################################################
#
#   LinearAlgebraicGroup, LinearAlgebraicGroupElem
#
###############################################################################

@doc raw"""
    LinearAlgebraicGroup

See [`linear_algebraic_group(::RootSystem)`](@ref) for the constructor.
"""
@attributes mutable struct LinearAlgebraicGroup <: Group
  #n::Int #GLn in which the LAG is embedded
  #type::Symbol # :SL, :SO, :Sp, :GL, :Other
  rs::RootSystem #in the future to be replaced by RootDatum
  G::MatrixGroup #the actual group
  k::Field #the field over which the group is defined

  function LinearAlgebraicGroup(R::RootSystem, G::MatrixGroup, k::Field)
    return new(R, G, k)
  end
end

function root_system(LAG::LinearAlgebraicGroup)
  return LAG.rs
end

function degree(LAG::LinearAlgebraicGroup)
  return degree(LAG.G)
end

@doc raw"""
    LinearAlgebraicGroupElem
"""

@attributes mutable struct LinearAlgebraicGroupElem <: GroupElem
  parent::LinearAlgebraicGroup
  mat::MatrixGroupElem #the actual eleemnt

  function LinearAlgebraicGroupElem(parent::LinearAlgebraicGroup, MGE::MatrixGroupElem)
    #add checks here
    return new(parent, MGE)
  end
end

function parent(a::LinearAlgebraicGroupElem)
  return a.parent
end
