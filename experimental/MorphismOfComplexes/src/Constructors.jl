function lift_morphism_through_free_resolutions(phi::ModuleFPHom;
    domain_resolution::ComplexOfMorphisms=free_resolution(domain(phi)).C,
    codomain_resolution::ComplexOfMorphisms=free_resolution(codomain(phi)).C
  )
  factory = LiftingMorphismsThroughResolution(phi)

  return MorphismOfComplexes(domain_resolution, codomain_resolution, factory, 0,
                             left_bound = -2, right_bound = -1, 
                             extends_left = false, extends_right = true)
end
