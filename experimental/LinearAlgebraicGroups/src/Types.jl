###############################################################################
#
#   LinearAlgebraicGroup, LinearAlgebraicGroupElem
#
###############################################################################

@doc raw"""
    LinearAlgebraicGroup

See [`linear_algebraic_group(::RootSystem)`](@ref) for the constructor.
"""
@attributes mutable struct LinearAlgebraicGroup<: Group
  #n::Int #GLn in which the LAG is embedded
  #type::Symbol # :SL, :SO, :Sp, :GL, :Other
  root_system::RootSystem #in the future to be replaced by RootDatum
  G::MatrixGroup #the actual group

  function LinearAlgebraicGroup(type::Symbol, n::Int)  
    if type == :A
      R = root_system(:A, n)
      G = special_linear_group(n + 1, QQ)
      LAG = new(R,G)
    else
      error("Only type A is implemented so far.")
    end
    return LAG
  end
end




@doc raw"""
    LinearAlgebraicGroupElem
"""
