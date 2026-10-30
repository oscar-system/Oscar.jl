########################################################################
# Direct sums of complexes
#
# This realizes direct sums of complexes of modules over a common ring 
# `R` in a lazy way. The anticipated use case are strands of complexes 
# which are composed from monomial bases for the 'fine grading'. There 
# a single strand for a coarse degree `alpha` consists of a whole 
# bunch of direct sums of strands for the fine grading for fine degrees 
# `d` with coarse degree `alpha`. Such direct sums are huge and usually 
# not all terms in a complex are needed. Moreover, it is crucial that 
# the functoriality is accessible in an economic way. 
########################################################################

### Production of the chains
struct DirectSumChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  R::Ring
  summands::Vector{<:AbsHyperComplex}

  function DirectSumChainFactory{CT}(R::Ring, summands::Vector{<:AbsHyperComplex}) where {CT}
    return new{CT}(R, summands)
  end
end

function (fac::DirectSumChainFactory)(self::AbsHyperComplex, i::Tuple)
  res = _direct_sum(fac.R, [can_compute_index(s, i) ? s[i] : _zero_module(fac.R) for s in fac.summands])
  return res
end

function can_compute(fac::DirectSumChainFactory, self::AbsHyperComplex, i::Tuple)
  return true
end

### Production of the morphisms 
struct DirectSumMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType} end

function (::DirectSumMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  fac = chain_factory(self)
  dom = self[i]
  I = collect(i)
  I[p] += direction(self, p) == :chain ? -1 : 1
  j = Tuple(I)
  cod = self[j]
  maps = OFPModuleHom{<:OFPModule, <:OFPModule, Nothing}[]
  for (k, s) in enumerate(fac.summands)
    if !can_compute_map(s, p, i)
      domk = domain(canonical_injection(dom, k))
      codk = domain(canonical_injection(cod, k))
      push!(maps, hom(domk, codk, elem_type(codk)[zero(cod_k) for _ in 1:ngens(domk)]))
      continue
    end
    push!(maps, map(s, p, i))
  end
  return _direct_sum(fac.R, maps; domain=dom, codomain=cod)
end

# We use the internal function `_direct_sum` instead of `direct_sum`. By default it 
# just deflects to `direct_sum`. But since we can not guarantee the latter to be streamlined 
# in its behavior for different types of modules and coefficient rings, this deviation allows 
# us to wrap any existing functionality so that it will adhere to the interface used here. 
# Note that, in particular, we need to pass the common ring `R` as an additional argument 
# to be able to catch the edge cases where zero modules/morphisms are provided. 

# Direct sums of morphisms. Keyword arguments allow to specify domain and codomain. 
function _direct_sum(R::Ring, phis::Vector{T}; 
    domain=_direct_sum(R, [domain(phi) for phi in phis]), 
    codomain=_direct_sum(R, [codomain(phi) for phi in phis])
  ) where {T<:OFPModuleHom}
  return direct_sum(phis; domain, codomain)
end

# By default we use `OFPModule`s here. In case no summands are provided, these internal 
# functions create the zero module according to the type of the ring. If you want this 
# to behave differently for your application, overwrite this for your type of rings. 
_zero_module(R::Ring) = FreeMod(R, 0)
_zero_module(R::MPolyDecRing) = graded_free_module(R, elem_type(grading_group(R))[])

function _direct_sum(R::Ring, summands::Vector{T}) where {T <: OFPModule}
  is_empty(summands) && return _zero_module(R)
  return direct_sum(summands; task=:none)
end

# In case the type of modules is not recognized by the input vector, throw an 
# error. The programmer should take care that their code is sufficiently type stable.
function _direct_sum(R::Ring, summands::Vector)
  @assert is_empty(summands) "non-empty list with no useful type detected"
  return _zero_module(R)
end

function _direct_sum(R::Ring, phis::Vector{<:OFPModuleHom{<:OFPModule, <:OFPModule, Nothing}}; 
    domain::OFPModule=_direct_sum(R, [domain(phi) for phi in phis]), 
    codomain::OFPModule=_direct_sum(R, [codomain(phi) for phi in phis])
  )
  img_gens = sizehint!(Vector{elem_type(codomain)}[], length(phis))
  for (k, phi) in enumerate(phis)
    inc_cod = canonical_injection(codomain, k)
    @assert Oscar.domain(inc_cod) === Oscar.codomain(phi)
    ig2 = images_of_generators(phi)
    ig3 = [inc_cod(v) for v in ig2]
    push!(img_gens, ig3)
  end
  return hom(domain, codomain, is_empty(img_gens) ? elem_type(codomain)[] : reduce(vcat, img_gens))
end

function can_compute(::DirectSumMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return true
end

### The concrete struct
@doc raw"""
    DirectSumComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
    
Direct sum of complexes of `OFPModule`s over a common ring `R`. In particular, this realizes
the induced morphisms. In order to be able to catch edge cases for terms with an empty list 
of summands, the ring itself needs to be provided to the constructor. 
"""
@attributes mutable struct DirectSumComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function DirectSumComplex(
      R::Ring, dirs::Vector{Symbol}, summands::Vector{<:AbsHyperComplex{CT, MT}}
    ) where {CT <: OFPModule, MT <: OFPModuleHom}
    return DirectSumComplex{CT, MT}(R, dirs, summands)
  end

  function DirectSumComplex{CT, MT}(
      R::Ring, dirs::Vector{Symbol}, summands::Vector
    ) where {CT <: OFPModule, MT <: OFPModuleHom}
    d = length(dirs)
    if !is_empty(summands)
      s = first(summands)::AbsHyperComplex
      @assert all(dim(ss) == d for ss in summands[2:end])
      @assert all(all(direction(ss, p) == dirs[p] for p in 1:dim(s)) for ss in summands[2:end])
    end
    chain_fac = DirectSumChainFactory{CT}(R, summands)
    map_fac = DirectSumMapFactory{MT}()

    internal_complex = HyperComplex(d, chain_fac, map_fac, dirs)
    return new{CT, MT}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::DirectSumComplex) = c.internal_complex

# dirty catch of the edge case with an empty list of summands 
# where the type parameter of the `summands` vector might not be specified
function DirectSumComplex(
    R::Ring, dirs::Vector{Symbol}, summands::Vector
  )
  @assert isempty(summands) "non-empty list of summands but no useful type of the entries of that list; try preparing your list in a type-stable way or specify the type parameters explicitly"
  ET = elem_type(R)
  return DirectSumComplex(R, dirs, AbsHyperComplex{OFPModule{ET}, OFPModuleHom}[])
end

# method for the user's convenience
function direct_sum(R::Ring, list::Vector{<:AbsHyperComplex})
  is_empty(list) && return DirectSumComplex{OFPModule{elem_type(R)}, OFPModuleHom}(R, Symbol[], list)
  dirs = Symbol[direction(first(list), i) for i in 1:dim(first(list))]
  return DirectSumComplex{OFPModule{elem_type(R)}, OFPModuleHom}(R, dirs, list)
end

