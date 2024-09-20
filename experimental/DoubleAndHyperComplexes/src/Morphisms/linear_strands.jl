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
  d::Int
  map_to_original::AbsHyperComplexMorphism

  function LinearStrandComplex(
      orig::AbsHyperComplex{ChainType, MorphismType};
      degree::Int = 0
    ) where {ChainType, MorphismType}
    chain_fac = LinearStrandChainFactory(orig, degree)
    map_fac = LinearStrandMapFactory(orig)
    # TODO: Do we want to also allow higher dimensional complexes here?
    @assert dim(orig) == 1 "complex must be 1-dimensional"
    @assert all(k->direction(orig, k) == :chain, 1:dim(orig)) "complex must be a chain complex in all directions"

    upper_bounds = Union{Nothing, Int}[(has_upper_bound(orig, i) ? upper_bound(orig, i) : nothing) for i in 1:dim(orig)]
    lower_bounds = Union{Nothing, Int}[(has_lower_bound(orig, i) ? lower_bound(orig, i) : nothing) for i in 1:dim(orig)]
    internal_complex = HyperComplex(dim(orig), chain_fac, map_fac, 
                                    [direction(orig, i) for i in 1:dim(orig)], 
                                    upper_bounds = upper_bounds,
                                    lower_bounds = lower_bounds
                                   )
    return new{ChainType, MorphismType}(internal_complex, orig, degree)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::LinearStrandComplex) = c.internal_complex
original_complex(c::LinearStrandComplex) = c.original_complex
inclusion_to_original_complex(c::LinearStrandComplex) = c.map_to_original
degree(c::LinearStrandComplex) = c.d

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

### User facing constructor
function linear_strand(C::AbsHyperComplex, degree::Int=0)
  result = LinearStrandComplex(C, degree=degree)
  result.map_to_original = LinearStrandInclusion(result)
  return result, result.map_to_original
end

function linear_strand(C::AbsHyperComplex, g::FinGenAbGroupElem)
  return linear_strand(C, g[1])
end

### Additional constructors for compatibility with what's already there:
# TODO: Adjust once we have decided which structures to use in the long run!
function linear_strand(C::ComplexOfMorphisms, degree::Int=0)
  return linear_strand(SimpleComplexWrapper(C), degree)
end

function linear_strand(C::FreeResolution, degree::Int=0)
  rng = range(C.C)
  # TODO: Using the operator [first(rng):-1:0] on C.C directly does not 
  # cut off the lower end. Is this a bug? 
  return linear_strand(SimpleComplexWrapper(C.C)[0:first(rng)], degree)
end

### Associated functionality for cokernels

### Production of the chains
struct LinearStrandComplementChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  str::LinearStrandComplex
  maps_from_original::Dict{<:Tuple, Any}

  function LinearStrandComplementChainFactory(
      str::LinearStrandComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    maps_from_original = Dict{Tuple, Any}()
    return new{ChainType}(str, maps_from_original)
  end
end

function (fac::LinearStrandComplementChainFactory)(self::AbsHyperComplex, i::Tuple)
  str = fac.str
  orig = original_complex(str)
  
  F_full = orig[i]::FreeMod
  @assert is_graded(F_full)
  S = base_ring(F_full)
  @assert is_standard_graded(S) "module and ring must be standard graded"
  G = grading_group(S)
  offset = zero(G) + sum(i)*G[1] + degree(str)*G[1]
  min_ind = [k for k in 1:rank(F_full) if degree(F_full[k]) == offset]
  comp = [k for k in 1:rank(F_full) if !(k in min_ind)]
  F = graded_free_module(S, elem_type(G)[degree(F_full[k]) for k in comp])
  map = hom(F_full, F, elem_type(F)[(k in comp ? F[findfirst(==(k), comp)] : zero(F)) for k in 1:rank(F_full)])
  fac.maps_from_original[i] = map
  return F
end

function can_compute(fac::LinearStrandComplementChainFactory, self::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.str, i)
end

### Production of the morphisms 
struct LinearStrandComplementMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  str::LinearStrandComplex

  function LinearStrandComplementMapFactory(
      str::LinearStrandComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    return new{MorphismType}(str)
  end
end

function (fac::LinearStrandComplementMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  str = fac.str
  orig = original_complex(str)
  orig_map = map(orig, p, i)
  dom = self[i]
  next = _codomain_index(self, p, i)
  cod = self[next]
  from_orig_dom = chain_factory(self).maps_from_original[i]
  from_orig_cod = chain_factory(self).maps_from_original[next]
  g = [preimage(from_orig_dom, x) for x in gens(dom)]
  img_gens = from_orig_cod.(orig_map.(g))
  return hom(dom, cod, img_gens)
end

function can_compute(fac::LinearStrandComplementMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(original_complex(fac.str), p, i)
end

### The concrete struct
@attributes mutable struct LinearStrandComplementComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  linear_strand::LinearStrandComplex{ChainType, MorphismType}

  function LinearStrandComplementComplex(
      str::LinearStrandComplex{ChainType, MorphismType}
    )  where {ChainType, MorphismType}
    chain_fac = LinearStrandComplementChainFactory(str)
    map_fac = LinearStrandComplementMapFactory(str)

    upper_bounds = Union{Int, Nothing}[(has_upper_bound(str, i) ? upper_bound(str, i) : nothing) for i in 1:dim(str)]
    lower_bounds = Union{Int, Nothing}[(has_lower_bound(str, i) ? lower_bound(str, i) : nothing) for i in 1:dim(str)]
    internal_complex = HyperComplex(dim(str), chain_fac, map_fac, [:chain for i in 1:dim(str)],
                                   upper_bounds = upper_bounds,
                                   lower_bounds = lower_bounds
                                  )
    return new{ChainType, MorphismType}(internal_complex, str)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::LinearStrandComplementComplex) = c.internal_complex
original_strand(c::LinearStrandComplementComplex) = c.linear_strand
original_complex(c::LinearStrandComplementComplex) = original_complex(original_strand(c))

### The projection map for the cokernel of the inclusion
struct LinearStrandQuotientProjectionFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  comp::LinearStrandComplementComplex

  function LinearStrandQuotientProjectionFactory(
      comp::LinearStrandComplementComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    return new{MorphismType}(comp)
  end
end

function (fac::LinearStrandQuotientProjectionFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  comp = fac.comp
  comp[i] # fill the cache
  return chain_factory(comp).maps_from_original[i]
end

function can_compute(fac::LinearStrandQuotientProjectionFactory, self::AbsHyperComplexMorphism, i::Tuple)
  can_compute_index(original_complex(fac.comp), i)
end


@attributes mutable struct LinearStrandQuotientProjection{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, LinearStrandQuotientProjection{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  comp::LinearStrandComplementComplex

  function LinearStrandQuotientProjection(
      comp::LinearStrandComplementComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    map_factory = LinearStrandQuotientProjectionFactory(comp)
    orig = original_complex(comp)

    internal_morphism = HyperComplexMorphism(orig, comp, map_factory, 
                                             cached=true, offset=[0 for i in 1:dim(orig)])
    return new{typeof(orig), typeof(comp), MorphismType}(internal_morphism, comp)
  end
end

underlying_morphism(phi::LinearStrandQuotientProjection) = phi.internal_morphism


### user facing constructor
function cokernel(inc::LinearStrandInclusion)
  str = domain(inc)
  orig = codomain(inc)
  d = degree(str)
  result = LinearStrandComplementComplex(str)
  pr = LinearStrandQuotientProjection(result)
  return result, pr
end
