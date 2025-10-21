#For generating the Database entries

@attributes mutable struct MatroidRealizations{BaseRingType, RingType} <: AbsAffineScheme{BaseRingType, RingType}
  matroid::Matroid 
  realization_space::MatroidRealizationSpace
  dim_r::Int
  selfproj_realization::MatroidRealizationSpaceSelfProjecting
  dim_s::Int
  equal::Bool
end
