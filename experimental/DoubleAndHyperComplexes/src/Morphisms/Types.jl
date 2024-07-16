abstract type AbsHyperComplexMorphism{DomainType<:AbsHyperComplex, CodomainType<:AbsHyperComplex, MorphismType, SelfType} <: Map{DomainType, CodomainType, SetMap, SelfType} end

# Map interface
domain(phi::AbsHyperComplexMorphism) = domain(underlying_morphism(phi))
codomain(phi::AbsHyperComplexMorphism) = codomain(underlying_morphism(phi))

# Get the morphism starting from D[i]
getindex(phi::AbsHyperComplexMorphism, i::Tuple) = underlying_morphism(phi)[i]
getindex(phi::AbsHyperComplexMorphism, i::Int...) = getindex(phi, i)

can_compute_index(phi::AbsHyperComplexMorphism, i::Tuple) = can_compute_index(underlying_morphism(phi), i)
can_compute_index(phi::AbsHyperComplexMorphism, i::Int...) = can_compute_index(phi, i)

has_index(phi::AbsHyperComplexMorphism, i::Tuple) = has_index(underlying_morphism(phi), i)
has_index(phi::AbsHyperComplexMorphism, i::Int) = has_index(phi, i)

# Get the offset k with the morphisms going from D[i] -> C[i+k].
offset(phi::AbsHyperComplexMorphism) = offset(underlying_morphism(phi))

### Factories for the production of morphisms 
abstract type HyperComplexMorphismFactory{MorphismType} end

function (fac::HyperComplexMorphismFactory)(Phi::AbsHyperComplexMorphism, i::Tuple)
  error("production of the $i-th entry not implemented for factories of type $(typeof(fac))")
end

function can_compute(fac::HyperComplexMorphismFactory, Phi::AbsHyperComplexMorphism, i::Tuple)
  error("deciding whether the $i-th entry can be produced by factories of type $(typeof(fac)) not implemented")
end


### Minimal concrete type for morphisms of hypercomplexes
mutable struct HyperComplexMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, HyperComplexMorphism{DomainType, CodomainType, MorphismType}} 
  domain::DomainType
  codomain::CodomainType
  fac::HyperComplexMorphismFactory{MorphismType}
  offset::Vector{Int}
  is_complete::Union{Bool, Nothing}
  map_cache::Dict{<:Tuple, <:MorphismType}


  function HyperComplexMorphism(
      dom::AbsHyperComplex, cod::AbsHyperComplex, 
      fac::HyperComplexMorphismFactory{MorphismType};
      cached::Bool=true, offset::Vector{Int}=[0 for i in 1:dim(dom)]
    ) where {MorphismType}
    if cached
      map_cache = Dict{Tuple, MorphismType}()
      return new{typeof(dom), typeof(cod), MorphismType}(dom, cod, fac, offset, nothing, map_cache)
    else
      return new{typeof(dom), typeof(cod), MorphismType}(dom, cod, fac, offset, nothing)
    end
  end
end

domain(phi::HyperComplexMorphism) = phi.domain
codomain(phi::HyperComplexMorphism) = phi.codomain

function getindex(phi::HyperComplexMorphism{<:Any, <:Any, MorphismType}, i::Tuple) where {MorphismType}
  # No caching deflects to production
  !isdefined(phi, :map_cache) && return phi.fac(phi, i)::MorphismType

  # The case of cached results
  haskey(phi.map_cache, i) && return phi.map_cache[i]
  can_compute(phi.fac, phi, i) || error("index out of bounds")
  result = phi.fac(phi, i)::MorphismType
  phi.map_cache[i] = result
  return result
end

offset(phi::HyperComplexMorphism) = phi.offset

function can_compute_index(phi::HyperComplexMorphism, i::Tuple)
  isdefined(phi, :map_cache) && haskey(phi.map_cache, i) && return true
  @assert length(i) == dim(domain(phi))
  !can_compute_index(domain(phi), i) && return false
  !can_compute_index(codomain(phi), Tuple(collect(i) + offset(phi))) && return false
  return can_compute(phi.fac, phi, i)
end

### Abstract types for morphisms of simple complexes

abstract type AbsSimpleComplexMorphism{DomainType<:AbsSimpleComplex, CodomainType<:AbsSimpleComplex, MorphismType, SelfType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, SelfType} end

### Widen the signature for the user's convenience
getindex(phi::AbsSimpleComplexMorphism, i::Int) = phi[(i,)]
can_compute_index(phi::AbsSimpleComplexMorphism, i::Int) = can_compute_index(phi, (i,))
has_index(phi::AbsSimpleComplexMorphism, i::Int) = has_index(phi, (i,))

mutable struct SimpleComplexMorphismWrapper{DomainType<:AbsSimpleComplex, CodomainType<:AbsSimpleComplex, MorphismType} <: AbsSimpleComplexMorphism{DomainType, CodomainType, MorphismType, SimpleComplexMorphismWrapper{DomainType, CodomainType, MorphismType}}
  phi::AbsHyperComplexMorphism

  function SimpleComplexMorphismWrapper(
      phi::AbsHyperComplexMorphism{D, C, M}
    ) where {D <:AbsSimpleComplex, C <:AbsSimpleComplex, M}
    return new{D, C, M}(phi)
  end
end

underlying_morphism(phi::SimpleComplexMorphismWrapper) = phi.phi

### Simplified complexes
@attributes mutable struct SimplifiedComplex{ChainType, MorphismType} <: AbsSimpleComplex{ChainType, MorphismType}
  complex::AbsHyperComplex{ChainType, MorphismType}
  original_complex::AbsHyperComplex{ChainType, MorphismType}
  map_to_original::AbsHyperComplexMorphism
  map_from_original::AbsHyperComplexMorphism
  equalizing_homotopy::AbsHyperComplexMorphism

  # The maps have to be set by the external constructor. 
  # They have this object as (co-)domain, so they can not be set 
  # here. 
  function SimplifiedComplex(
      complex::AbsHyperComplex{ChainType, MorphismType},
      original_complex::AbsHyperComplex{ChainType, MorphismType}
      ) where {ChainType, MorphismType}
    @assert dim(original_complex) == dim(complex) == 1 "complexes must be of dimension one"
    return new{ChainType, MorphismType}(complex, original_complex)
  end
end

underlying_complex(c::SimplifiedComplex) = c.complex
map_to_original_complex(c::SimplifiedComplex) = c.map_to_original
map_from_original_complex(c::SimplifiedComplex) = c.map_from_original
original_complex(c::SimplifiedComplex) = c.original_complex

@attributes mutable struct SimpleFreeResolution{ChainType, MorphismType, ModuleType} <: AbsSimpleComplex{ChainType, MorphismType}
  M::ModuleType # The original module to be resolved
  underlying_complex::AbsHyperComplex{ChainType, MorphismType}
  augmentation_map::AbsHyperComplexMorphism # A map to the 1-term complex  0 ← M ← 0 in degree 0

  # The augmentation map has this object as its domain, so it needs 
  # to be set by the external constructor.
  function SimpleFreeResolution(
      M::T,
      complex::AbsHyperComplex{ChainType, MorphismType},
    ) where {ChainType, MorphismType, T}
    return new{ChainType, MorphismType, T}(M, complex)
  end
end

underlying_complex(c::SimpleFreeResolution) = c.underlying_complex
original_module(c::SimpleFreeResolution) = c.M

