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
  summands::Vector{<:AbsHyperComplex}

  function DirectSumChainFactory{CT}(summands::Vector{<:AbsHyperComplex}) where {CT}
    return new{CT}(summands)
  end
end

function (fac::DirectSumChainFactory)(self::AbsHyperComplex, i::Tuple)
  ind = findfirst(can_compute_index(a, i) for a in fac.summands)
  isnothing(ind) && error("can not compute this index")
  sample = fac.summands[ind][i]
  res, _ = direct_sum([can_compute_index(s, i) ? s[i] : _zero_object(sample) for s in fac.summands])
  return res
end

function can_compute(fac::DirectSumChainFactory, self::AbsHyperComplex, i::Tuple)
  return any(can_compute_index(a, i) for a in fac.summands)
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
  return direct_sum(maps; domain=dom, codomain=cod)
end

# We need to be able to create dummy zero maps for those cases 
# where one of the summands does not have objects in the respective 
# homological degrees and hence can not produce maps themselves.
# Since there is no canonical way to do this, we request the user 
# to fill this gap here according to the type they wish to use.
function _zero_map(domain, codomain)
  error("not implemented for arguments of type $(typeof(domain)); please overwrite this internal for the direct sums of complexes, see source code for details")
end

# `DirectSumComplex` produces an entry at index `i` if and only if 
# there is at least one summand which can produce an entry at index `i`. 
# For all others a zero object needs to be created to fill the spots 
# in the structure of the direct sum. Since there is no standard way 
# to do this, the user must overwrite this function for their types. 
function _zero_object(sample)
  error("not implemented for arguments of type $(typeof(sample)); please overwrite this internal for the direct sums of complexes, see source code for details")
end

# For instance, for `OFPModule`s we do the following: 
function _zero_object(M::OFPModule)
  is_graded(M) && return graded_free_module(base_ring(M), grading_group(M)[])
  return FreeMod(base_ring(M), 0)
end

function _zero_map(domain::OFPModule, codomain::OFPModule)
  @assert base_ring(domain) === base_ring(codomain) "incompatible rings in domain and codomain"
  return hom(domain, codomain, elem_type(codomain)[zero(codomain) for _ in 1:ngens(domain)])
end

function can_compute(::DirectSumMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  return any(can_compute_map(a, p, i) for a in chain_factory(self).summands)
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

  function DirectSumComplex{CT, MT}(
      summands::Vector
    ) where {CT, MT}
    d = is_empty(summands) ? 0 : dim(first(summands))
    dirs = is_empty(summands) ? Symbol[] : Symbol[direction(first(summands), i) for i in 1:d]
    if !is_empty(summands)
      s = first(summands)::AbsHyperComplex
      @assert all(dim(ss) == d for ss in summands[2:end])
      @assert all(all(direction(ss, p) == dirs[p] for p in 1:dim(s)) for ss in summands[2:end])
    end
    chain_fac = DirectSumChainFactory{CT}(summands)
    map_fac = DirectSumMapFactory{MT}()

    internal_complex = HyperComplex(d, chain_fac, map_fac, dirs)
    return new{CT, MT}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::DirectSumComplex) = c.internal_complex

# method for the user's convenience
function direct_sum(list::Vector{<:AbsHyperComplex})
  error("no reasonable type for objects and morphisms could be inferred from the argument list; try making the type of the list more informative, or call the internal constructor directly")
end

function direct_sum(list::Vector{<:AbsHyperComplex{CT}}) where {CT<:OFPModule}
  return DirectSumComplex{CT, OFPModuleHom}(list)
end

