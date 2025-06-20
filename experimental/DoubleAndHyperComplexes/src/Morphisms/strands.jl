########################################################################
# Strands of hyper complexes of graded free modules.
#
# Suppose S = 𝕜[x₀,…,xₙ] is a standard graded ring and C a 
# hypercomplex of modules over it with all morphisms of degree 
# zero. Then on every module M in C there is a degree-0-part 
# and the boundary maps restrict to these to form a hypercomplex 
# of 𝕜-modules. 
#
# We produce these strand complexes, together with their inclusion 
# maps into their original complex of S-modules.
########################################################################

### Production of the chains
struct StrandChainFactory{ChainType<:ModuleFP} <: HyperComplexChainFactory{ChainType}
  orig::AbsHyperComplex
  d::Union{Int, FinGenAbGroupElem}

  function StrandChainFactory(
      orig::AbsHyperComplex{ChainType}, d::Union{Int, FinGenAbGroupElem}
    ) where {ChainType<:ModuleFP}
    return new{FreeMod}(orig, d) # TODO: Specify the chain type better
  end
end

### Production of the morphisms 
struct StrandMorphismFactory{MorphismType<:ModuleFPHom} <: HyperComplexMapFactory{MorphismType}
  orig::AbsHyperComplex
  d::Union{Int, FinGenAbGroupElem}
  monomial_mappings::Dict{<:Tuple{<:Tuple, Int}, <:Map}

  function StrandMorphismFactory(orig::AbsHyperComplex, d::Union{Int, FinGenAbGroupElem})
    monomial_mappings = Dict{Tuple{Tuple, Int}, Map}()
    return new{FreeModuleHom}(orig, d, monomial_mappings)
  end
end

### The concrete struct
@attributes mutable struct StrandComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  original_complex::AbsHyperComplex
  d::Union{Int, FinGenAbGroupElem}
  inclusion_map::AbsHyperComplexMorphism

  function StrandComplex(
      orig::AbsHyperComplex{ChainType, MorphismType}, d::Union{Int, FinGenAbGroupElem}
    ) where {ChainType <: ModuleFP, MorphismType <: ModuleFPHom}
    chain_fac = StrandChainFactory(orig, d)
    map_fac = StrandMorphismFactory(orig, d)

    internal_complex = HyperComplex(dim(orig), 
                                    chain_fac, map_fac, Symbol[direction(orig, i) for i in 1:dim(orig)],
                                    upper_bounds = [(has_upper_bound(orig, p) ? upper_bound(orig, p) : nothing) for p in 1:dim(orig)],
                                    lower_bounds = [(has_lower_bound(orig, p) ? lower_bound(orig, p) : nothing) for p in 1:dim(orig)],
                                   )
    # Assuming that ChainType and MorphismType are provided by the input
    return new{FreeMod, FreeModuleHom}(internal_complex, orig, d)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::StrandComplex) = c.internal_complex

degree(c::StrandComplex) = c.d
original_complex(c::StrandComplex) = c.original_complex
inclusion_map(c::StrandComplex) = c.inclusion_map

### The morphism for the strand's inclusion
struct StrandInclusionMorphismFactory{MorphismType<:ModuleFPHom} <: HyperComplexMorphismFactory{MorphismType}
  strand::AbsHyperComplex

  function StrandInclusionMorphismFactory(strand::AbsHyperComplex)
    return new{FreeModuleHom}(strand)
  end
end

function (fac::StrandInclusionMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  orig = codomain(self)
  strand = domain(self)
  dom = strand[i]
  cod = orig[i]
  R = base_ring(cod)
  to = hom(dom, cod, elem_type(cod)[prod(x^k for (x, k) in zip(gens(R), e); init=one(R))*cod[i] for (e, i) in all_exponents(cod, degree(strand))])
  return to
end

function can_compute(fac::StrandInclusionMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(domain(self), i)
end


@attributes mutable struct StrandInclusionMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, StrandInclusionMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function StrandInclusionMorphism(strand::StrandComplex)
    map_factory = StrandInclusionMorphismFactory(strand)
    dom = strand
    cod = original_complex(strand)

    # Assuming that the domain `dom` and the codomain `cod` have 
    # been extracted from the input
    internal_morphism = HyperComplexMorphism(dom, cod, map_factory)
    # Assuming that `MorphismType` has been extracted from the input
    return new{typeof(dom), typeof(cod), FreeModuleHom}(internal_morphism)
  end
end

underlying_morphism(phi::StrandInclusionMorphism) = phi.internal_morphism

