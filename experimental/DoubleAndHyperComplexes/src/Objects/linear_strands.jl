########################################################################
# Linear strands according to Eisenbud: The Geometry of Syzygies, 
# Chapter 7.
########################################################################

### Production of the chains
struct LinearStrandChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  orig::AbsHyperComplex
  d::Int # the degree of this linear strand
  maps_to_orig::Dict{<:Tuple, Any}  # for building the maps later.

  function LinearStrandChainFactory(
      orig::AbsHyperComplex{ChainType},
      d::Int
    ) where {ChainType}
    maps_to_orig = Dict{Tuple, Any}()
    return new{ChainType}(orig, d, maps_to_orig)
  end
end

function (fac::LinearStrandChainFactory)(self::AbsHyperComplex, i::Tuple)
  F_full = fac.orig[i]::FreeMod
  @assert is_graded(F_full)
  S = base_ring(F_full)
  @assert is_standard_graded(S) "module and ring must be standard graded"
  G = grading_group(S)
  offset = zero(G) + sum(i)*G[1] + fac.d*G[1]
  min_ind = [k for k in 1:rank(F_full) if degree(F_full[k]) == offset]
  F = graded_free_module(S, [offset for i in 1:length(min_ind)])
  map = hom(F, F_full, elem_type(F_full)[F_full[i] for i in min_ind])
  fac.maps_to_orig[i] = map
  return F
end

function can_compute(fac::LinearStrandChainFactory, self::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.orig, i)
end

### Production of the morphisms 
struct LinearStrandMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  orig::AbsHyperComplex

  function LinearStrandMapFactory(
      orig::AbsHyperComplex{ChainType, MorphismType}
    ) where {ChainType, MorphismType}
    return new{MorphismType}(orig)
  end
end

function (fac::LinearStrandMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  orig_map = map(fac.orig, p, i)
  dom = self[i]
  next = _codomain_index(self, p, i)
  cod = self[next]
  to_orig_dom = chain_factory(self).maps_to_orig[i]
  to_orig_cod = chain_factory(self).maps_to_orig[next]
  img_gens = [preimage(to_orig_cod, orig_map(to_orig_dom(v))) for v in gens(dom)]
  return hom(dom, cod, img_gens)
end

function _codomain_index(self::AbsHyperComplex, p::Int, i::Tuple)
  I = collect(i)
  I[p] = (direction(self, p) == :chain ? I[p] - 1 : I[p] + 1)
  return Tuple(I)
end


function can_compute(fac::LinearStrandMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(fac.orig, p, i)
end

### The concrete struct has been postponed to ../Morphisms/linear_strands.jl
# due to its dependency on various types of morphisms which are only introduced 
# later.

