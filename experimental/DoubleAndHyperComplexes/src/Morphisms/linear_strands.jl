########################################################################
# Linear strands according to Eisenbud: The Geometry of Syzygies, 
# Chapter 7.
########################################################################

struct LinearStrandInclusionFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  function LinearStrandInclusionFactory(::Type{MorphismType}) where {MorphismType}
    return new{MorphismType}()
  end
end

function (fac::LinearStrandInclusionFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  str = domain(self)
  orig = codomain(self)
  str[i] # fill the cache
  return chain_factory(str).maps_to_orig[i]
end

function can_compute(fac::LinearStrandInclusionFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(domain(self), i)
end

# Postponed concrete type for linear strands (need AbsHyperComplexMorphism)
@attributes mutable struct LinearStrandComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  original_complex::AbsHyperComplex{ChainType, MorphismType}
  map_to_original::AbsHyperComplexMorphism

  function LinearStrandComplex(
      orig::AbsHyperComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    chain_fac = LinearStrandChainFactory(orig)
    map_fac = LinearStrandMapFactory(orig)

    upper_bounds = Union{Nothing, Int}[(has_upper_bound(orig, i) ? upper_bound(orig, i) : nothing) for i in 1:dim(orig)]
    lower_bounds = Union{Nothing, Int}[(has_lower_bound(orig, i) ? lower_bound(orig, i) : nothing) for i in 1:dim(orig)]
    internal_complex = HyperComplex(dim(orig), chain_fac, map_fac, 
                                    [direction(orig, i) for i in 1:dim(orig)], 
                                    upper_bounds = upper_bounds,
                                    lower_bounds = lower_bounds
                                   )
    return new{ChainType, MorphismType}(internal_complex, orig)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::LinearStrandComplex) = c.internal_complex
original_complex(c::LinearStrandComplex) = c.original_complex
inclusion_to_original_complex(c::LinearStrandComplex) = c.map_to_original

### Concrete type for the inclusion
@attributes mutable struct LinearStrandInclusion{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, LinearStrandInclusion{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function LinearStrandInclusion(
      str::LinearStrandComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    map_factory = LinearStrandInclusionFactory(MorphismType)

    internal_morphism = HyperComplexMorphism(str, original_complex(str), 
                                             map_factory, cached=true, 
                                             offset=[0 for i in 1:dim(str)])
    return new{typeof(str), typeof(original_complex(str)), MorphismType}(internal_morphism)
  end
end

underlying_morphism(phi::LinearStrandInclusion) = phi.internal_morphism

# User facing constructor
function linear_strand(C::AbsHyperComplex)
  result = LinearStrandComplex(C)
  result.map_to_original = LinearStrandInclusion(result)
  return result, result.map_to_original
end

