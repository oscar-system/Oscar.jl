include("realization_space.jl")
include("SelfProjectingMatroids.jl")
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
export is_selfprojecting
export selfprojecting_realization_space
