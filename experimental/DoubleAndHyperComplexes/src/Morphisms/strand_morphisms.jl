#= 
# Induced maps on strands of complexes
#
# Say one has a morphism of complexes of graded modules
# `ϕ : C* → D*` and a degree `e`. One wishes to compute 
# the induced morphism 
# 
#   `ϕₑ : Cₑ* → Dₑ*`
#
# in the degree-`e`-part.
=#

struct InducedStrandMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  phi::AbsHyperComplexMorphism
  d::FinGenAbGroupElem
  domain::StrandComplex
  codomain::StrandComplex

  function InducedStrandMorphismFactory(
      phi::AbsHyperComplexMorphism,
      d::FinGenAbGroupElem,
      domain::StrandComplex,
      codomain::StrandComplex
    )
    return new{ModuleFPHom}(phi, d, domain, codomain)
  end
end

function (fac::InducedStrandMorphismFactory)(self::AbsHyperComplexMorphism, I::Tuple)
  II = collect(I)
  JJ = II + offset(fac.phi)
  J = Tuple(JJ)
  dom = fac.domain[I]
  cod = fac.codomain[J]
  orig_dom = domain(fac.phi)[I]
  orig_cod = codomain(fac.phi)[J]
  orig_map = fac.phi[I]

  # Use a dictionary for fast mapping of the monomials to the 
  # generators of `cod`.
  cod_dict = Dict{Tuple{Vector{Int}, Int}, elem_type(cod)}(m=>cod[k] for (k, m) in enumerate(all_exponents(orig_cod, fac.d)))
  # Hashing of FreeModElem's can not be assumed to be non-trivial. Hence we use the exponents directly.
  img_gens_res = elem_type(cod)[]
  R = base_ring(orig_dom)
  vv = gens(R)
  for (e, i) in all_exponents(orig_dom, fac.d) # iterate through the generators of `dom`
    m = prod(x^k for (x, k) in zip(vv, e); init=one(R))*orig_dom[i]
    v = orig_map(m) # map the monomial
    # take preimage of the result using the previously built dictionary.
    # TODO: Iteration over the terms of v is VERY slow due to its suboptimal implementation.
    # We have to iterate manually. This saves us roughly 2/3 of the memory consumption and 
    # it also runs three times as fast. 
    w = zero(cod)
    for (i, b) in coordinates(v)
      #g = orig_cod[i]
      w += sum(c*cod_dict[(n, i)] for (c, n) in zip(AbstractAlgebra.coefficients(b), AbstractAlgebra.exponent_vectors(b)); init=zero(cod))
    end
    push!(img_gens_res, w)
  end
  return hom(dom, cod, img_gens_res)
end

function can_compute(fac::InducedStrandMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.domain, i) && can_compute_index(fac.codomain, Tuple(collect(i) + offset(fac.phi)))
end


@attributes mutable struct InducedStrandMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, InducedStrandMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function InducedStrandMorphism(
      phi::AbsHyperComplexMorphism,
      d::FinGenAbGroupElem;
      domain::StrandComplex = strand(Oscar.domain(phi), d)[1],
      codomain::StrandComplex = strand(Oscar.codomain(phi), d)[1]
    )
    @assert original_complex(domain) === Oscar.domain(phi)
    @assert original_complex(codomain) === Oscar.codomain(phi)
    @assert degree(domain) == degree(codomain) == d

    map_factory = InducedStrandMorphismFactory(phi, d, domain, codomain)
    
    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, 
                                             cached=true, offset=offset(phi))
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::InducedStrandMorphism) = phi.internal_morphism

function strand(
    phi::AbsHyperComplexMorphism,
    d::FinGenAbGroupElem;
    domain::StrandComplex = strand(Oscar.domain(phi), d)[1],
    codomain::StrandComplex = strand(Oscar.codomain(phi), d)[1]
  )
  return InducedStrandMorphism(phi, d; domain, codomain)
end

degree(c::InducedStrandMorphism) = degree(domain(c))

