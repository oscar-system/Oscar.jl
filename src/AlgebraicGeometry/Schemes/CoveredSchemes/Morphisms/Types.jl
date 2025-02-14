
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

### Compatibility
underlying_morphism(f::CoveredSchemeMorphism) = f

function CoveredSchemeMorphism(f::CoveredSchemeMorphism)
  return f
end



##############################################################################
## Concrete Type for normalization
## very similar to CoveredSchemeMorphism, but allowing disjoint handling
## of disjoint components
##############################################################################
@doc raw"""
    NormalizationMorphism{
                  DomainType<:AbsCoveredScheme,
                  CodomainType<:AbsCoveredScheme
      } <:AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 Nothing,
                                 NormalizationMorphism,
                                }
A datastructure to encode normalizations of covered schemes.

It is described as the morphism from the new scheme to the original one, containing
information on the decomposition of the new scheme into disjoint components.
(This is the type of the return value of  `normalization(X::AbsCoveredScheme)`.)
"""
@attributes mutable struct NormalizationMorphism{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme
   } <:AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 Nothing,
                                 NormalizationMorphism,
                                }

  underlying_morphism::CoveredSchemeMorphism
  inclusions::Vector{<:AbsCoveredSchemeMorphism}

  function NormalizationMorphism(
      f::CoveredSchemeMorphism,
      inclusions::Vector{<:AbsCoveredSchemeMorphism};
      check::Bool=true
    ) 
    @check is_normal(X) "not a normalization morphism"
    @assert all(inc->codomain(inc) === domain(f), inclusions) "domains and codomains do not match"
    ret_value = new{typeof(domain(f)),typeof(codomain(f))}(f,inclusions)
    return ret_value    
  end

  function NormalizationMorphism(
      f::CoveredSchemeMorphism;
      check::Bool=true
    )
    @check is_normal(X) "not a normalization morphism"
    ret_value = new{typeof(domain(f)),typeof(codomain(f))}(f,[identity_map(X)])
    return ret_value    
  end

end
