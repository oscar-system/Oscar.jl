function underlying_morphism(inc::PrincipalOpenEmbedding)
  return inc.inc
end

function complement_equations(inc::PrincipalOpenEmbedding)
  return inc.complement_equations
end

function image(inc::PrincipalOpenEmbedding)
  isdefined(inc, :image) && return inc.image::PrincipalOpenSubset

  result = PrincipalOpenSubset(codomain(inc), complement_equations(inc))
  inc.image = result
  return result
end

function inverse_on_image(inc::PrincipalOpenEmbedding)
  if isdefined(inc, :inverse_on_image) 
    return inc.inverse_on_image::AffineSchemeMor
  end

  X = domain(inc)
  Y = codomain(inc)
  U = image(inc)

  if ambient_coordinate_ring(X) === ambient_coordinate_ring(U)
    result = morphism(U, X, gens(OO(X)); check=false)
    inc.inverse_on_image = result
    return result
  end

  res = morphism(X, U, pullback(inc).(gens(OO(Y))), check=false)
  result = inverse(res)
  inc.inverse_on_image = result
  return result
end

