#For generating the Database entries

struct MatroidRealizations
  name::String #name will be of the form "rk_3_n_6_index_1"
  matroid::Matroid 
  rank::Int
  length_groundset::Int
  realization_space::MatroidRealizationSpace
  dim_r::Int
  selfprojecting_realization_space::Union{MatroidRealizationSpaceSelfProjecting,Nothing} # in case that the computation of the selfprojecting realization space did not terminate the value will be nothing, as for dim_s and equal
  dim_s::Union{Int,Nothing}
  equality_of_realizationspaces::Union{Bool,Nothing}
end

matroid(MR::MatroidRealizations) = MR.matroid
realization_space(MR::MatroidRealizations) = MR.realization_space
selfprojecting_realization_space(MR::MatroidRealizations) = MR.selfprojecting_realization_space
dim_r(MR::MatroidRealizations) = MR.dim_r
dim_s(MR::MatroidRealizations) = MR.dim_s
equality_of_realizationspaces(MR::MatroidRealizations) = MR.equality_of_realizationspaces
name(MR::MatroidRealizations) = MR.name
length_groundset(MR::MatroidRealizations) = MR.length_groundset
rank(MR::MatroidRealizations) = MR.rank

function Base.show(io::IO, ::MIME"text/plain", MR::MatroidRealizations)
  io = Oscar.pretty(io)
  print(io, "The matroid is of rank ", MR.rank, " on ", MR.length_groundset, " elements.\n")
  show(io, MIME("text/plain"),realization_space(MR))
  print(io, "\n")
  if !isnothing(selfprojecting_realization_space(MR))
    show(io, MIME("text/plain"),selfprojecting_realization_space(MR))
    print(io, "\nThe closures of the realization space and the self-projecting realization space are ")
    if equality_of_realizationspaces(MR)
      print(io, "equal.")
    elseif !equality_of_realizationspaces(MR)
      print(io, "not equal.")
    end
  else
    print(io, "The computation of the self-projecting realization space did not terminate.")
  end
end

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