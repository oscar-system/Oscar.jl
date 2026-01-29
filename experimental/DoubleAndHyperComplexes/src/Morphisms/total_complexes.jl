#= 
# Induced morphisms on total complexes
#
# Say one has a morphism of suitably bounded complexes `ϕ : C* → D*`. 
# One wishes to compute the induced morphism on the total complexes 
#
#   tot(ϕ) : tot(C*) → tot(D*).
#
=#

struct TotalComplexMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  phi::AbsHyperComplexMorphism
  dom::TotalComplex
  cod::TotalComplex

  function TotalComplexMorphismFactory(
      phi::AbsHyperComplexMorphism,
      dom::TotalComplex,
      cod::TotalComplex
    )
    return new{ModuleFPHom}(phi, dom, cod)
  end
end

function (fac::TotalComplexMorphismFactory)(self::AbsHyperComplexMorphism, I::Tuple)
  # fill the cache
  dom = fac.dom[I]
  J = Tuple(first(I) + first(offset(self)))
  cod = fac.cod[J]
  dom_inds = indices_in_summand(fac.dom, first(I)) 
  cod_inds = indices_in_summand(fac.dom, first(J))
  result = hom(dom, cod, elem_type(cod)[zero(cod) for _ in 1:ngens(dom)])
  # We only need to treat the "overlapping" part of the index ranges.
  for i in dom_inds
    j = Tuple(collect(i) + offset(fac.phi))
    j in cod_inds || continue
    result += compose(projection(fac.dom, i), compose(fac.phi[i], injection(fac.cod, j)))
  end
  return result
end

function can_compute(fac::TotalComplexMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.dom, i) && can_compute_index(fac.cod, i)
end


@attributes mutable struct TotalComplexMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, TotalComplexMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function TotalComplexMorphism(
      phi::AbsHyperComplexMorphism;
      domain::TotalComplex = total_complex(Oscar.domain(phi)),
      codomain::TotalComplex = total_complex(Oscar.codomain(phi))
    )
    @assert original_complex(domain) === Oscar.domain(phi)
    @assert original_complex(codomain) === Oscar.codomain(phi)
    map_factory = TotalComplexMorphismFactory(phi, domain, codomain)
    
    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, cached=true, offset=[sum(offset(phi); init=0)])
    # Assuming that the types have been extracted from the input
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::TotalComplexMorphism) = phi.internal_morphism

function total_complex(
    phi::AbsHyperComplexMorphism;
    domain::TotalComplex = total_complex(Oscar.domain(phi)),
    codomain::TotalComplex = total_complex(Oscar.codomain(phi))
  )
  return TotalComplexMorphism(phi; domain, codomain)
end

