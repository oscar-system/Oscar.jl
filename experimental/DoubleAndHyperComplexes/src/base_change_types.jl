# change of base rings for complexes of morphisms

### Production of the chains
struct BaseChangeChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  phi::Any # Whatever serves as a base change map
  orig::AbsHyperComplex
  red_map_cache::Dict{Tuple, Any}

  function BaseChangeChainFactory(phi::Any, orig::AbsHyperComplex)
    # TODO: Can we be more specific about the type?
    return new{ModuleFP}(phi, orig, Dict{Tuple, Any}())
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
    return new{ModuleFPHom}(phi, orig)
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
  internal_complex::HyperComplex{ChainType, MorphismType}

  function BaseChangeComplex(phi::Any, orig::AbsHyperComplex{CT, MT}) where {CT<:ModuleFP, MT<:ModuleFPHom}
    chain_fac = BaseChangeChainFactory(phi, orig)
    map_fac = BaseChangeMapFactory(phi, orig)

    # Assuming d is the dimension of the new complex
    d = dim(orig)
    internal_complex = HyperComplex(d, chain_fac, map_fac, [direction(orig, i) for i in 1:d])
    # Assuming that ChainType and MorphismType are provided by the input
    return new{ModuleFP, ModuleFPHom}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::BaseChangeComplex) = c.internal_complex

