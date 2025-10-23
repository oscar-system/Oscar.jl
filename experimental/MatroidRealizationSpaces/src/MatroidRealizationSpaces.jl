include("realization_space.jl")
include("SelfProjectingMatroids.jl")
include("for_database.jl")
include("serialization.jl")

### Methods for realization spaces of matroids
function underlying_scheme(RS::MatroidRealizationSpace{BRT, RT}) where {BRT<:Ring, RT<:MPolyQuoLocRing}
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme::AffineScheme{BRT, RT}

  P = ambient_ring(RS)::MPolyRing
  I = defining_ideal(RS)::MPolyIdeal
  U = MPolyPowersOfElement(P, P.(inequations(RS)))::MPolyPowersOfElement
  RS.underlying_scheme = spec(P, I, U)
  return RS.underlying_scheme::AffineScheme{BRT, RT}
end

# Exports
export inequations
export is_realizable
export realization
export realization_matrix
export realization_space
export defining_ideal

#Exports from SelfProjectingMatroids.jl
export is_selfprojecting
export satisfies_disjointbasisproperty
export selfprojecting_realization_space
export MatroidRealizationSpace_SelfProj
export underlying_scheme
export selfproj_realization_ideal
export selfproj_realization_matrix
export MatroidRealizationSpace_SelfProj
export inequations
export defining_ideal
export underlying_scheme
export selfproj_realization_ideal
export dimension
export basisminors

#Exports from for_database.jl
export matroid 
export realization_space
export selfproj_realization_space
export dim_r
export dim_s
export equal
export name
export length_groundset
export rk
export MatroidRealizations
