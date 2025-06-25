#= 
# Given a morphism of 1-dimensional complexes `ϕ : C* → D*` one 
# wishes to compute `ϕ : C[a:b]* → D[a':b']*`
# where `a:b` is a range and `a':b'` the same range shifted by the `offset` of `ϕ`.
#
# Note: At the moment we only support 1-dimensional complexes 
# and ordinary ranges a:b with positive orientation.
=#

struct ViewMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  phi::AbsHyperComplexMorphism
  rng::UnitRange{Int}
  domain::HyperComplexView
  codomain::HyperComplexView

  function ViewMorphismFactory(
      phi::AbsHyperComplexMorphism,
      rng::UnitRange{Int},
      domain::HyperComplexView,
      codomain::HyperComplexView
    )
    return new{ModuleFPHom}(phi, rng, domain, codomain)
  end
end

function (fac::ViewMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  return fac.phi[i]
end

function can_compute(fac::ViewMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return first(i) in fac.rng
end


@attributes mutable struct ViewMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, ViewMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function ViewMorphism(
      phi::AbsHyperComplexMorphism,
      rng::UnitRange{Int};
      domain::HyperComplexView = Oscar.domain(phi)[rng],
      codomain::HyperComplexView = Oscar.codomain(phi)[first(rng)+first(offset(phi)):last(rng)+first(offset(phi))]
    )
    @assert Oscar.domain(phi) === original_complex(domain)
    @assert Oscar.codomain(phi) === original_complex(codomain)
    # TODO: Check compatibility of ranges
    map_factory = ViewMorphismFactory(phi, rng, domain, codomain)

    internal_morphism = HyperComplexMorphism(domain, codomain, map_factory, cached=false, offset=offset(phi))
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::ViewMorphism) = phi.internal_morphism

function getindex(
    phi::AbsHyperComplexMorphism, rng::UnitRange{Int};
    domain::HyperComplexView = Oscar.domain(phi)[rng],
    codomain::HyperComplexView = Oscar.codomain(phi)[first(rng)+first(offset(phi)):last(rng)+first(offset(phi))]
  )
  return ViewMorphism(phi, rng; domain, codomain)
end
