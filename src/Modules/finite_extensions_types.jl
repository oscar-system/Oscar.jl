mutable struct FiniteExtension{DomainType, CodomainType} <: Map{DomainType, CodomainType, Map, FiniteExtension}
  phi::MPolyAnyMap{DomainType, CodomainType, Nothing}
  fiber_ideal::Ideal
  fiber_quotient::MPolyQuoRing
  basis::Vector
  codomain_as_module::SubquoModule
  interp

  function FiniteExtension(
      phi::MPolyAnyMap{DomainType, CodomainType, Nothing}
    ) where {DomainType, CodomainType}
    return new{DomainType, CodomainType}(phi)
  end
end

