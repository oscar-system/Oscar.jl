# Given a 1-dimensional hypercomplex `C` of graded modules with 
# homogeneous (co-)boundary maps, this constructs the same complex,
# but with all modules shifted so that the morphisms become homogeneous 
# of degree zero. 
#
# Note that to this end the complex needs to be connected in the sense 
# that there must not be zero modules or maps sitting between any non-zero 
# modules. 
#
# Furthermore, one module needs to be given from which the levelling to 
# degree zero starts. By default, this is the codomain of the last 
# nonzero (co-)boundary map (depending on the direction of the complex).

### Production of the chains
struct DegreeZeroChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  orig::AbsHyperComplex{ChainType}
  deltas::Dict{<:Tuple, Any}

  function DegreeZeroChainFactory(orig::AbsHyperComplex{ChainType}) where {ChainType}
    return new{ChainType}(orig, Dict{Tuple, Any}())
  end
end

function (fac::DegreeZeroChainFactory)(self::AbsHyperComplex, i::Tuple)
  dom = fac.orig[i]

  if !can_compute_map(self, 1, i)
    fac.deltas[i] = zero(grading_group(dom))
    return dom
  end

  if iszero(dom) 
    G = grading_group(dom)
    fac.deltas[i] = zero(G)
    return dom
  end

  is_graded(dom) || error("module must be graded")
  S = base_ring(dom)
  G = grading_group(S)
  j = _codomain_index(fac.orig, 1, i)
  cod = self[j] # Make sure the cache for the deltas is filled
  phi = map(fac.orig, 1, i)

  delta = degree(phi)
  d = fac.deltas[j] - delta
  fac.deltas[i] = d
  # 
  # delta = elem_type(G)[]
  # for g in gens(dom)
  #   d = degree(g)
  #   h = phi(g)
  #   iszero(h) && error("generator must not map to zero")
  #   e = degree(phi(g))
  #   push!(delta, d - e)
  # end
  #
  # @assert all(k->degree(k) == first(delta), delta) "degree shifts must coincide"
  # 
  result = twist(dom, d)
  return result
end

function can_compute(fac::DegreeZeroChainFactory, self::AbsHyperComplex, i::Tuple)
  if direction(fac.orig, 1) == :chain && has_lower_bound(fac.orig, 1)
    return can_compute_index(fac.orig, i) && first(i) >= lower_bound(fac.orig, 1)
  else
    return can_compute_index(fac.orig, i) && first(i) <= upper_bound(fac.orig, 1)
  end
end

### Production of the morphisms 
struct DegreeZeroMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  orig::AbsHyperComplex

  function DegreeZeroMapFactory(orig::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}
    return new{MorphismType}(orig)
  end
end

function (fac::DegreeZeroMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  dom = self[i]
  orig_dom = fac.orig[i]
  j = _codomain_index(self, p, i)
  cod = self[j]
  cod_orig = fac.orig[j]
  phi_orig = map(fac.orig, p, i)
  img_gens_orig = images_of_generators(phi_orig)
  new_imgs = [cod(coordinates(g)) for g in img_gens_orig]
  return hom(dom, cod, new_imgs, check=false)
end

function can_compute(fac::DegreeZeroMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  can_compute_map(fac.orig, p, i) || return false
  phi = map(fac.orig, p, i)
  is_graded(phi) || return false
  any(iszero, images_of_generators(phi)) && return false
  return true
end

### The concrete struct
@attributes mutable struct DegreeZeroComplex{ChainType, MorphismType} <: AbsSimpleComplex{ChainType, MorphismType} 
  internal_complex::AbsHyperComplex{ChainType, MorphismType}
  orig::AbsHyperComplex{ChainType, MorphismType}

  function DegreeZeroComplex(orig::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}
    dim(orig) == 1 || error("complex must have dimension one")
    if direction(orig, 1) == :chain 
      has_lower_bound(orig, 1) || error("chain complex must be bounded from below")
    else
      has_upper_bound(orig, 1) || error("cochain complex must be bounded from above")
    end
    chain_fac = DegreeZeroChainFactory(orig)
    map_fac = DegreeZeroMapFactory(orig)

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(
        1, chain_fac, map_fac, [direction(orig, 1)], 
        upper_bounds = Union{Int, Nothing}[has_upper_bound(orig, 1) ? upper_bound(orig) : Nothing],
        lower_bounds = Union{Int, Nothing}[has_lower_bound(orig, 1) ? lower_bound(orig) : Nothing]
       )
    # Assuming that ChainType and MorphismType are provided by the input
    return new{ChainType, MorphismType}(internal_complex, orig)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::DegreeZeroComplex) = c.internal_complex

