function getindex(moc::MorphismOfComplexes{<:Any, <:Any, MorphismType}, i::Int) where {MorphismType}
  has_left_bound(moc) && !extends_left(moc) && (i >= left_bound(moc) || error("index out of range"))
  has_right_bound(moc) && !extends_right(moc) && (i <= right_bound(moc) || error("index out of range"))

  haskey(moc.morphisms, i) && return moc.morphisms[i]::MorphismType
  result = moc.morphism_factory(moc, i)::MorphismType
  moc.morphisms[i] = result
  return result
end

function (fac::LiftingMorphismsThroughResolution)(phi::AbsMorphismOfComplexes, i::Int)
  i == -1 && return fac.original_map
  i == -2 && return hom(domain(phi)[-2], codomain(phi)[-2], elem_type(codomain(phi)[-2])[]) # the zero morphism

  # Build up the maps recursively via liftings
  psi = phi[i-1]
  dom_map = map(domain(phi), i)
  cod_map = map(codomain(phi), i)
  x = gens(domain(dom_map))
  y = psi.(dom_map.(x))
  z = [preimage(cod_map, u) for u in y]
  return hom(domain(dom_map), domain(cod_map), z)
end

