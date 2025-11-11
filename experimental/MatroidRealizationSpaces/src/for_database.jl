#For generating the Database entries

struct MatroidRealizations
  name::String #name will be of the form "rk_3_n_6_index_1"
  matroid::Matroid 
  rk::Int
  length_groundset::Int
  realization_space::MatroidRealizationSpace
  dim_r::Int
  selfproj_realization_space::MatroidRealizationSpaceSelfProjecting
  dim_s::Int
  equal::Bool
end

matroid(MR::MatroidRealizations) = MR.matroid
realization_space(MR::MatroidRealizations) = MR.realization_space
selfproj_realization_space(MR::MatroidRealizations) = MR.selfproj_realization_space
dim_r(MR::MatroidRealizations) = MR.dim_r
dim_s(MR::MatroidRealizations) = MR.dim_s
equal(MR::MatroidRealizations) = MR.equal
name(MR::MatroidRealizations) = MR.name
length_groundset(MR::MatroidRealizations) = MR.length_groundset
rk(MR::MatroidRealizations) = MR.rk

# #Export
# export matroid 
# export realization_space
# export selfproj_realization
# export dim_r
# export dim_s
# export equal
# export name
# export length_groundset
# export rk