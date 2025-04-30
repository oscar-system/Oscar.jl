# change of base rings for complexes of morphisms

### Production of the chains
struct BaseChangeChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  phi::Any # Whatever serves as a base change map
  orig::AbsHyperComplex
  red_map_cache::Dict{Tuple, Any}

  function BaseChangeChainFactory(phi::Any, orig::AbsHyperComplex)
    # TODO: Can we be more specific about the type?
    return new{SparseFPModule}(phi, orig, Dict{Tuple, Any}())
  end
end

function (fac::BaseChangeChainFactory)(self::AbsHyperComplex, i::Tuple)
  res, m = change_base_ring(fac.phi, fac.orig[i])
  fac.red_map_cache[i] = m
  return res
end

function can_compute(fac::BaseChangeChainFactory, self::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.orig, i)
end

### Production of the morphisms 
struct BaseChangeMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  phi::Any # Whatever serves as a base change map
  orig::AbsHyperComplex

  function BaseChangeMapFactory(phi::Any, orig::AbsHyperComplex)
    return new{SparseFPModuleHom}(phi, orig)
  end
end

function (fac::BaseChangeMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  f = map(fac.orig, p, i)
  next = _codomain_index(self, p, i)
  # Fill the cache
  self[i]
  self[next]
  dom_bc = chain_factory(self).red_map_cache[i]
  cod_bc = chain_factory(self).red_map_cache[next]
  res = change_base_ring(fac.phi, f; domain_base_change=dom_bc, codomain_base_change=cod_bc)
end

function can_compute(fac::BaseChangeMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(fac.orig, p, i)
end

### The concrete struct
@attributes mutable struct BaseChangeComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  orig::AbsHyperComplex
  internal_complex::HyperComplex{ChainType, MorphismType}
  red_map::AbsHyperComplexMorphism

  function BaseChangeComplex(phi::Any, orig::AbsHyperComplex{CT, MT}) where {CT<:SparseFPModule, MT<:SparseFPModuleHom}
    chain_fac = BaseChangeChainFactory(phi, orig)
    map_fac = BaseChangeMapFactory(phi, orig)

    d = dim(orig)
    upper_bounds = Vector{Union{Int, Nothing}}([(has_upper_bound(orig, i) ? upper_bound(orig, i) : nothing) for i in 1:d])
    lower_bounds = Vector{Union{Int, Nothing}}([(has_lower_bound(orig, i) ? lower_bound(orig, i) : nothing) for i in 1:d])
    internal_complex = HyperComplex(d, chain_fac, map_fac, 
                                    [direction(orig, i) for i in 1:d],
                                    lower_bounds = lower_bounds,
                                    upper_bounds = upper_bounds
                                   )
    return new{SparseFPModule, SparseFPModuleHom}(orig, internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::BaseChangeComplex) = c.internal_complex

########################################################################
# Type for the base change morphism
########################################################################
struct BaseChangeMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  cod::BaseChangeComplex

  function BaseChangeMorphismFactory(comp::BaseChangeComplex{CT, MT}) where {CT, MT}
    # TODO: Can we do more about the type?
    return new{SparseFPModuleHom}(comp)
  end
end

function (fac::BaseChangeMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  # fill the cache
  fac.cod[i]

  return chain_factory(fac.cod).red_map_cache[i]
end

function can_compute(fac::BaseChangeMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.cod, i)
end


@attributes mutable struct BaseChangeMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, BaseChangeMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function BaseChangeMorphism(cod::BaseChangeComplex)
    map_factory = BaseChangeMorphismFactory(cod)
    dom = cod.orig

    internal_morphism = HyperComplexMorphism(dom, cod, map_factory, cached=true, offset=[0 for i in 1:dim(dom)])
    return new{typeof(cod.orig), typeof(cod), SparseFPModuleHom}(internal_morphism)
  end
end

underlying_morphism(phi::BaseChangeMorphism) = phi.internal_morphism

