### Production of the chains
struct DirectSumChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  R::Ring
  summands::Vector{<:AbsHyperComplex}

  function DirectSumChainFactory(R::Ring, summands::Vector{<:AbsHyperComplex{CT}}) where {CT}
    return new{CT}(R, summands)
  end
end

function (fac::DirectSumChainFactory)(self::AbsHyperComplex, i::Tuple)
  res = _direct_sum(fac.R, [s[i] for s in fac.summands])
  return res
end

function can_compute(fac::DirectSumChainFactory, self::AbsHyperComplex, i::Tuple)
  is_empty(fac.summands) && return true
  return can_compute_index(first(fac.summands), i)
  return all(can_compute_index(s, i) for s in fac.summands)
end

### Production of the morphisms 
struct DirectSumMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function DirectSumMapFactory(::Type{MorphismType}) where MorphismType
    return new{MorphismType}()
  end
end

function (::DirectSumMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  fac = chain_factory(self)
  dom = self[i]
  I = collect(i)
  I[p] += direction(self, p) == :chain ? -1 : 1
  j = Tuple(I)
  cod = self[j]
  return _direct_sum(fac.R, [map(s, p, i) for s in fac.summands]; domain=dom, codomain=cod)
end

function _direct_sum(R::Ring, phis::Vector{T}; 
    domain=_direct_sum(R, [domain(phi) for phi in phis]), 
    codomain=_direct_sum(R, [codomain(phi) for phi in phis])
  ) where {T<:OFPModuleHom}
  return direct_sum(phis; domain, codomain)
end

_zero_module(R::Ring) = FreeMod(R, 0)
_zero_module(R::MPolyDecRing) = graded_free_module(R, elem_type(grading_group(R))[])

function _direct_sum(R::Ring, summands::Vector{T}) where {T <: OFPModule}
  is_empty(summands) && return _zero_module(R)
  return direct_sum(summands; task=:none)
end

function _direct_sum(R::Ring, summands::Vector)
  @assert is_empty(summands) "non-empty list with no useful type detected"
  return _zero_module(R)
end

function direct_sum(phis::Vector{<:OFPModuleHom{<:OFPModule, <:OFPModule, Nothing}}; 
    domain::OFPModule=direct_sum([domain(phi) for phi in phis])[1],
    codomain::OFPModule=direct_sum([codomain(phi) for phi in phis])[1]
  )
  img_gens = elem_type(codomain)[]
  for (k, phi) in enumerate(phis)
    pr_dom = canonical_projection(domain, k)
    @assert Oscar.codomain(pr_dom) === Oscar.domain(phi)
    inc_cod = canonical_injection(codomain, k)
    @assert Oscar.domain(inc_cod) === Oscar.codomain(phi)
    ig2 = images_of_generators(phi)
    ig3 = [inc_cod(v) for v in ig2]
    img_gens = vcat(img_gens, ig3)
  end
  return hom(domain, codomain, img_gens)
end

function _direct_sum(R::Ring, phis::Vector{<:OFPModuleHom{<:OFPModule, <:OFPModule, Nothing}}; 
    domain::OFPModule=_direct_sum(R, [domain(phi) for phi in phis]), 
    codomain::OFPModule=_direct_sum(R, [codomain(phi) for phi in phis])
  )
  img_gens = sizehint!(Vector{elem_type(codomain)}[], length(phis))
  for (k, phi) in enumerate(phis)
    #pr_dom = canonical_projection(domain, k)
    #@assert Oscar.codomain(pr_dom) === Oscar.domain(phi)
    inc_cod = canonical_injection(codomain, k)
    @assert Oscar.domain(inc_cod) === Oscar.codomain(phi)
    ig2 = images_of_generators(phi)
    ig3 = [inc_cod(v) for v in ig2]
    push!(img_gens, ig3)
  end
  return hom(domain, codomain, is_empty(img_gens) ? elem_type(codomain)[] : reduce(vcat, img_gens))
end

function can_compute(::DirectSumMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  fac = chain_factory(self)
  is_empty(fac.summands) && return true
  return can_compute_map(first(fac.summands), p, i)
  return all(can_compute_map(s, p, i) for s in fac.summands)
end

### The concrete struct
@attributes mutable struct DirectSumComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function DirectSumComplex(
      R::Ring, dirs::Vector{Symbol}, summands::Vector{<:AbsHyperComplex{CT, MT}}
    ) where {CT <: OFPModule, MT <: OFPModuleHom}
    d = length(dirs)
    if !is_empty(summands)
      s = first(summands)
      @assert all(dim(ss) == d for ss in summands[2:end])
      @assert all(all(direction(ss, p) == dirs[p] for p in 1:dim(s)) for ss in summands[2:end])
    end
    chain_fac = DirectSumChainFactory(R, summands)
    map_fac = DirectSumMapFactory(MT)

    internal_complex = HyperComplex(d, chain_fac, map_fac, dirs)
    return new{CT, MT}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::DirectSumComplex) = c.internal_complex

