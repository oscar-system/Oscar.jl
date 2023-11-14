include("realization_space.jl")

### Methods for realization spaces of matroids
function underlying_scheme(RS::MatroidRealizationSpace{BRT, RT}) where {BRT<:Ring, RT<:MPolyQuoLocRing}
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme::Spec{BRT, RT}

  P = ambient_ring(RS)::MPolyRing
  I = defining_ideal(RS)::MPolyIdeal
  U = MPolyPowersOfElement(P, P.(inequations(RS)))::MPolyPowersOfElement
  RS.underlying_scheme = Spec(P, I, U)
  return RS.underlying_scheme::Spec{BRT, RT}
end

