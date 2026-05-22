struct InducedBaseChangeMorphismFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType}
  bc::Any
  orig::AbsHyperComplexMorphism
  domain::BaseChangeComplex
  codomain::BaseChangeComplex
  run_checks::Bool
#=
  function InducedBaseChangeMorphismFactory{MorphismType}(bc::Any, 
                                                          orig::AbsHyperComplexMorphism, 
                                                          domain::BaseChangeComplex, 
                                                          codomain::BaseChangeComplex,
                                                          run_checks::Bool
                                                         )
    return new{MorphismType}(bc, orig, domain, codomain, run_checks)
  end
  =#
end

function (fac::InducedBaseChangeMorphismFactory)(self::AbsHyperComplexMorphism, i::Tuple)
  bc_dom = base_change_map(fac.domain)
  bc_cod = base_change_map(fac.codomain)
  if fac.run_checks
    @assert all(bc_dom[i](g) == gen(fac.domain[i], j) for (j, g) in enumerate(gens(domain(fac.orig[i])))) "generators do not map to generators on the domain base change map"
    @assert all(bc_cod[i](g) == gen(fac.codomain[i], j) for (j, g) in enumerate(gens(codomain(fac.orig[i])))) "generators do not map to generators on the codomain base change map"
  end
  img_gens = elem_type(fac.codomain[i])[bc_cod[i](v) for v in images_of_generators(fac.orig[i])]
  return hom(fac.domain[i], fac.codomain[i], img_gens)
end

function can_compute(fac::InducedBaseChangeMorphismFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.orig, i)
end


@attributes mutable struct InducedBaseChangeMorphism{DomainType, CodomainType, MorphismType} <: AbsHyperComplexMorphism{DomainType, CodomainType, MorphismType, InducedBaseChangeMorphism{DomainType, CodomainType, MorphismType}}
  internal_morphism::HyperComplexMorphism{DomainType, CodomainType, MorphismType}

  function InducedBaseChangeMorphism(
      bc::Any, orig::AbsHyperComplexMorphism;
      domain::AbsHyperComplex=change_base_ring(bc, domain(orig))[1], 
      codomain::AbsHyperComplex=change_base_ring(bc, codomain(orig))[1],
      check::Bool=true
    )
    map_factory = InducedBaseChangeMorphismFactory{ModuleFPHom}(bc, orig, domain, codomain, check)

    internal_morphism = HyperComplexMorphism(domain, codomain, 
                                             map_factory, cached=true, 
                                             offset=[0 for i in 1:dim(domain)])
    return new{typeof(domain), typeof(codomain), ModuleFPHom}(internal_morphism)
  end
end

underlying_morphism(phi::InducedBaseChangeMorphism) = phi.internal_morphism

