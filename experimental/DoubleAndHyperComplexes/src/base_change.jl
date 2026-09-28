function change_base_ring(phi::Any, C::AbsHyperComplex)
  res = BaseChangeComplex(phi, C)
  red_map = BaseChangeMorphism(res)
  res.red_map = red_map
  return res, red_map
end

function change_base_ring(
    bc::Any, phi::AbsHyperComplexMorphism;
    domain::BaseChangeComplex=change_base_ring(bc, domain(phi))[1],
    codomain::BaseChangeComplex=change_base_ring(bc, codomain(phi))[1],
    check::Bool=true
  )
  result = InducedBaseChangeMorphism(bc, phi; domain, codomain, check)
  return result, base_change_map(domain), base_change_map(codomain)
end

function base_change_map(comp::BaseChangeComplex)
  # The field is not set when using the internal constructor of the type.
  if !isdefined(comp, :red_map)
    comp.red_map = BaseChangeMorphism(comp)
  end
  return comp.red_map::BaseChangeMorphism
end

function original_complex(comp::BaseChangeComplex)
  return comp.orig
end

