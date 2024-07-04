
########################################################################
# Morphisms of covered schemes                                         #
########################################################################
abstract type AbsCoveredSchemeMorphism{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme,
    BaseMorphismType,
    CoveredSchemeMorphismType
   } <: SchemeMor{DomainType, CodomainType, CoveredSchemeMorphismType, BaseMorphismType}
end

########################################################################
# Concrete minimal type for morphisms of covered schemes               #
########################################################################
@attributes mutable struct CoveredSchemeMorphism{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme,
    BaseMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 CoveredSchemeMorphism,
                                 BaseMorphismType
                                }
  X::DomainType
  Y::CodomainType
  f::CoveringMorphism

  function CoveredSchemeMorphism(
      X::DomainType,
      Y::CodomainType,
      f::CoveringMorphism{<:Any, <:Any, BaseMorType};
      check::Bool=true
    ) where {
             DomainType<:AbsCoveredScheme,
             CodomainType<:AbsCoveredScheme,
             BaseMorType
            }
    @check is_refinement(domain(f), default_covering(X)) "covering not found in domain"
    @check is_refinement(codomain(f), default_covering(Y)) "covering not found in codomain"
    return new{DomainType, CodomainType, BaseMorType}(X, Y, f)
  end
end

