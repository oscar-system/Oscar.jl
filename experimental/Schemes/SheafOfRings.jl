@attribute mutable struct SheafOfRings{
    CoveredSchemeType<:CoveredScheme,
    CoveringType<:Covering,
    SpecType<:Spec,
    RingType<:AbstractAlgebra.Ring,
    RingHomType
  }
  X::CoveredSchemeType # the underlying scheme
  C::Covering # the original covering on which this sheaf had been defined
  rings::Dict{SpecType, RingType}
  identifications::Dict{Tuple{SpecType, SpecType}, RingHomType}
end

