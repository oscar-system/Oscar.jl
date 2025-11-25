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

function underlying_scheme(RS::MatroidRealizationSpace{BRT, RT}) where {BRT<:Ring, RT<:MPolyQuoRing}
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme::AffineScheme{BRT, RT}

  RS.underlying_scheme = spec(ambient_ring(RS)) # for some reason this is an MPolyQuoRing already
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
export selfprojecting_realization_ideal
export selfprojecting_realization_matrix
export inequations
export defining_ideal
export basisminors

#Exports from for_database.jl
export matroid 
export realization_space
export selfprojecting_realization_space
export dim_r
export dim_s
export equality_of_realizationspaces
export name
export length_groundset
export rank
export MatroidRealizations