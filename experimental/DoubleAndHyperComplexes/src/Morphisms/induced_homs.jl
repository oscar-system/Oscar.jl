#= 
# Induced complexes on `HomComplex`es.
#
# Say, one has morphisms of hypercomplexes `ϕ : C* → D*` and `ψ : E* → F`.
# Then we would like to be able to get all induced maps in the following 
# diagram:
#
#                   ϕ^*
#        Hom(D*, E*) → Hom(C*, E*)
#         ψ_* ↓            ↓ ψ_*
#        Hom(D*, F*) → Hom(C*, F*)
#                   ϕ^*
#        
# This is achieved by the co- and the contravariant induced homomorphisms.
=#
struct InducedContravariantMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  phi::AbsHyperComplexMorphism
  dom::HomComplex
  cod::HomComplex
  common_cod::AbsHyperComplex

  # Fields needed for production
  function InducedContravariantMorphismFactory(
      phi::AbsHyperComplexMorphism,
      D::AbsHyperComplex,
      dom::HomComplex,
      cod::HomComplex
    )
    return new{ModuleFPHom}(phi, dom, cod, D)
  end
end

function hom(phi::ModuleFPHom, M::ModuleFP; 
    domain::ModuleFP=hom(Oscar.codomain(phi), M),
    codomain::ModuleFP=hom(Oscar.domain(phi), M)
  )
  return lift_homomorphism_contravariant(domain, codomain, phi)
end

function (fac::InducedContravariantMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  dom = fac.dom
  cod = fac.cod
  common_cod = fac.common_cod
  n = dim(common_cod)
  m = dim(dom) - n
  a = i[1:m]
  b = i[m+1:m+n]
  na = Tuple([-i for i in a])
  return hom(fac.phi[na], common_cod[b]; domain=dom[i], codomain=cod[i])
end

function can_compute(fac::InducedContravariantMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  # Decide whether the outgoing map from index i can be computed
  dom = fac.dom
  cod = fac.cod
  common_cod = fac.common_cod
  n = dim(common_cod)
  m = dim(dom) - n
  a = i[1:m]
  b = i[m+1:m+n]
  na = Tuple([-i for i in a])
  return can_compute_index(fac.phi, na) && can_compute_index(common_cod, b)
end


@attributes mutable struct InducedContravariantMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, InducedContravariantMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  # Further specialized fields
  # ...

  function InducedContravariantMorphism(
      phi::AbsHyperComplexMorphism,
      D::AbsHyperComplex;
      domain::HomComplex = hom(Oscar.codomain(phi), D),
      codomain::HomComplex = hom(Oscar.domain(phi), D)
  )
    @assert Oscar.codomain(domain) === Oscar.codomain(codomain) === D
    @assert Oscar.domain(domain) === Oscar.codomain(phi)
    @assert Oscar.domain(codomain) === Oscar.domain(phi)
    map_factory = InducedContravariantMorphismFactory(phi, D, domain, codomain)

    # Assuming that the domain `dom` and the codomain `cod` have 
    # been extracted from the input
    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, cached=true, offset=[0 for i in 1:dim(domain)])
    # Assuming that the types have been extracted from the input
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::InducedContravariantMorphism) = phi.internal_morphism
  
function hom(
    phi::AbsHyperComplexMorphism, 
    D::AbsHyperComplex;
    domain::HomComplex = hom(Oscar.codomain(phi), D),
    codomain::HomComplex = hom(Oscar.domain(phi), D)
  )
  return InducedContravariantMorphism(phi, D; domain, codomain)
end

# The same for covariant induced maps
struct InducedCovariantMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  phi::AbsHyperComplexMorphism
  dom::HomComplex
  cod::HomComplex
  common_dom::AbsHyperComplex

  # Fields needed for production
  function InducedCovariantMorphismFactory(
      phi::AbsHyperComplexMorphism,
      D::AbsHyperComplex,
      dom::HomComplex,
      cod::HomComplex
    )
    return new{ModuleFPHom}(phi, dom, cod, D)
  end
end

function hom(M::ModuleFP, phi::ModuleFPHom; 
    domain::ModuleFP=hom(M, Oscar.domain(phi)),
    codomain::ModuleFP=hom(M, Oscar.codomain(phi))
  )
  return lift_homomorphism_covariant(domain, codomain, phi)
end

function (fac::InducedCovariantMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  dom = fac.dom
  cod = fac.cod
  common_dom = fac.common_dom
  n = dim(common_dom)
  m = dim(dom) - n
  a = i[1:n]
  b = i[n+1:m+n]
  na = Tuple([-i for i in a])
  return hom(common_dom[na], fac.phi[b]; domain=dom[i], codomain=cod[i])
end

function can_compute(fac::InducedCovariantMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  # Decide whether the outgoing map from index i can be computed
  dom = fac.dom
  cod = fac.cod
  common_dom = fac.common_dom
  n = dim(common_dom)
  m = dim(dom) - n
  a = i[1:n]
  b = i[n+1:m+n]
  na = Tuple([-i for i in a])
  return can_compute_index(fac.phi, b) && can_compute_index(common_dom, na)
end


@attributes mutable struct InducedCovariantMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, InducedCovariantMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}
  # Further specialized fields
  # ...

  function InducedCovariantMorphism(
      D::AbsHyperComplex,
      phi::AbsHyperComplexMorphism;
      domain::HomComplex = hom(D, Oscar.codomain(phi)),
      codomain::HomComplex = hom(D, Oscar.domain(phi))
  )
    @assert Oscar.domain(domain) === Oscar.domain(codomain) === D
    @assert Oscar.codomain(domain) === Oscar.domain(phi)
    @assert Oscar.codomain(codomain) === Oscar.codomain(phi)
    map_factory = InducedCovariantMorphismFactory(phi, D, domain, codomain)

    # Assuming that the domain `dom` and the codomain `cod` have 
    # been extracted from the input
    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, cached=true, offset=[0 for i in 1:dim(domain)])
    # Assuming that the types have been extracted from the input
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::InducedCovariantMorphism) = phi.internal_morphism
  
function hom(
    D::AbsHyperComplex,
    phi::AbsHyperComplexMorphism;
    domain::HomComplex = hom(D, Oscar.domain(phi)),
    codomain::HomComplex = hom(D, Oscar.codomain(phi))
  )
  return InducedCovariantMorphism(D, phi; domain, codomain)
end

