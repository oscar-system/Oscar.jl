### Production of the chains
struct BaseChangeChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  bc::Any
  orig::AbsHyperComplex
  bc_map_cache::Dict{Tuple, Map}

  function BaseChangeChainFactory{ChainType}(bc::Any, orig::AbsHyperComplex)
    return new{ChainType}(bc, orig, Dict{Tuple, Map}())
  end
end

function (fac::BaseChangeChainFactory)(self::AbsHyperComplex, i::Tuple)
  res, bc = change_base_ring(fac.bc, fac.orig[i])
  fac.bc_map_cache[i] = bc
  return res
end

function can_compute(fac::BaseChangeChainFactory, self::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.orig, i)
end

### Production of the morphisms 
struct BaseChangeMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  bc::Any
  orig::AbsHyperComplex

  function BaseChangeMapFactory{MorphismType}(bc::Any, orig::AbsHyperComplex)
    return new{MorphismType}(bc, orig)
  end
end

function (fac::BaseChangeMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  I = collect(i)
  I[p] += direction(self, p) == :chain ? -1 : 1
  j = Tuple(I)
  return change_base_ring(fac.bc, map(fac.orig, p, i); domain=self[i], codomain=self[j])
end

function can_compute(fac::BaseChangeMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(fac.orig, p, i)
end

### The concrete struct
@attributes mutable struct BaseChangeComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  bc::Any
  orig::AbsHyperComplex
  internal_complex::HyperComplex{ChainType, MorphismType}
  bc_map::AbsHyperComplexMorphism

  function BaseChangeComplex{ChainType, MorphismType}(bc::Any, orig::AbsHyperComplex)
    chain_fac = BaseChangeChainFactory{ChainType}(bc, orig)
    map_fac = BaseChangeMapFactory{MorphismType}(bc, orig)

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(d, chain_fac, map_fac, [direction(orig, i) for i in 1:d])
    # Assuming that ChainType and MorphismType are provided by the input
    return new{ChainType, MorphismType}(bc, orig, internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::BaseChangeComplex) = c.internal_complex

