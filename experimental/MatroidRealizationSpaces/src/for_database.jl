#For generating the Database entries

struct MatroidRealizations{BaseRingType, RingType} <: AbsAffineScheme{BaseRingType, RingType}
  name::String
  matroid::Matroid 
  rk::Int
  n::Int
  realization_space::MatroidRealizationSpace
  dim_r::Int
  selfproj_realization::MatroidRealizationSpaceSelfProjecting
  dim_s::Int
  equal::Bool
end
