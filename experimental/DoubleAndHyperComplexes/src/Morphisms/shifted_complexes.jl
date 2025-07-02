#= 
# Induced maps in shifted complexes
#
# Say one has a morphism `ϕ : C* → D*` and a shift `k`. Then one 
# would like to compute `ϕ[k] : C[k]* → D[k]*.
=#

struct ShiftedComplexMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  phi::AbsHyperComplexMorphism
  dom::ShiftedHyperComplex
  cod::ShiftedHyperComplex
  shift::Vector{Int}


  function ShiftedComplexMorphismFactory(
      phi::AbsHyperComplexMorphism,
      dom::ShiftedHyperComplex,
      cod::ShiftedHyperComplex
    )
    return new{ModuleFPHom}(phi, dom, cod, shift(dom))
  end
end

function (fac::ShiftedComplexMorphismFactory)(self::AbsHyperComplexMorphism, I::Tuple)
  j = Tuple([i for i in I] + fac.shift)
  return fac.phi[j]
end

function can_compute(fac::ShiftedComplexMorphismFactory, self::AbsHyperComplexMorphism, I::Tuple)
  j = Tuple([i for i in I] + fac.shift)
  return can_compute_index(fac.phi, j)
end


@attributes mutable struct ShiftedComplexMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, ShiftedComplexMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function ShiftedComplexMorphism(
      phi::AbsHyperComplexMorphism, d::Vector{Int};
      domain::ShiftedHyperComplex = shift(Oscar.domain(phi), d...),
      codomain::ShiftedHyperComplex = shift(Oscar.codomain(phi), d...)
    )
    @assert original_complex(domain) === Oscar.domain(phi)
    @assert original_complex(codomain) === Oscar.codomain(phi)
    @assert shift(domain) == shift(codomain) == d

    map_factory = ShiftedComplexMorphismFactory(phi, domain, codomain)

    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, cached=false, offset=offset(phi))
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::ShiftedComplexMorphism) = phi.internal_morphism

function shift(
    phi::AbsHyperComplexMorphism, d::Vector{Int};
    domain::ShiftedHyperComplex = shift(Oscar.domain(phi), d...),
    codomain::ShiftedHyperComplex = shift(Oscar.codomain(phi), d...)
  )
  return ShiftedComplexMorphism(phi, d; domain, codomain)
end

