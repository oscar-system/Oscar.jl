#For generating the Database entries

struct SelfProjectingMatroidRealizations
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

@doc raw"""
    matroid(MR::SelfProjectingMatroidRealizations)

Return the matroid of `MR`.
"""
matroid(MR::SelfProjectingMatroidRealizations) = MR.matroid

@doc raw"""
    realization_space(MR::SelfProjectingMatroidRealizations)

Return the realization space of `MR`.
"""
realization_space(MR::SelfProjectingMatroidRealizations) = MR.realization_space

@doc raw"""
    selfprojecting_realization_space(MR::SelfProjectingMatroidRealizations)

Return the self-projecting realization space of `MR`.
"""
selfprojecting_realization_space(MR::SelfProjectingMatroidRealizations) = MR.selfprojecting_realization_space

@doc raw"""
    dim_r(MR::SelfProjectingMatroidRealizations)

Return the (affine) dimension of the realization space of `MR`.
"""
dim_r(MR::SelfProjectingMatroidRealizations) = MR.dim_r

@doc raw"""
    dim_s(MR::SelfProjectingMatroidRealizations)

Return the (affine) dimension of the self-projecting realization space of `MR`."""
dim_s(MR::SelfProjectingMatroidRealizations) = MR.dim_s

@doc raw"""
    equality_of_realizationspaces(MR::SelfProjectingMatroidRealizations)

Return a boolean which states whether the self-projecting realization space  and the realization space of `MR` are equal."""
equality_of_realizationspaces(MR::SelfProjectingMatroidRealizations) = MR.equality_of_realizationspaces

@doc raw"""
    name(MR::SelfProjectingMatroidRealizations)

Return the  identifier of `MR` in the database."""
name(MR::SelfProjectingMatroidRealizations) = MR.name

@doc raw"""
    length_groundset(MR::SelfProjectingMatroidRealizations)

Return the size of the groundset of the matroid underlying `MR`."""
length_groundset(MR::SelfProjectingMatroidRealizations) = MR.length_groundset

@doc raw"""
    rank(MR::SelfProjectingMatroidRealizations)

Return the rank of the matroid underlying `MR`."""
rank(MR::SelfProjectingMatroidRealizations) = MR.rank

#function to deal with a nice display of the Matroid realizations
function Base.show(io::IO, ::MIME"text/plain", MR::SelfProjectingMatroidRealizations)
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
