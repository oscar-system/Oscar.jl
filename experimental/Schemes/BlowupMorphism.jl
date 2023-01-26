mutable struct BlowupMorphism{
                              CodomainType<:AbsCoveredScheme
                             }
  projective_bundle::AbsCoveredProjectiveScheme 
  codomain::CodomainType   # i.A. ein CoveredScheme
  center::IdealSheaf      # on codomain
  projection::AbsCoveredSchemeMor
  domain::AbsCoveredScheme # i.A. ein CoveredScheme
  exceptional_divisor::WeilDivisor

  function BlowupMorphism(
      IP::AbsCoveredProjectiveScheme,
      I::IdealSheaf
    )
    X = base_scheme(IP)
    X === scheme(I) || error("ideal sheaf not compatible with blown up variety")
    return new{typeof(X)}(IP, X, I)
  end
end

function domain(p::BlowupMorphism)
  if !isdefined(p, :domain)
    p.domain = covered_scheme(p.projective_bundle)
  end
  return p.domain
end

codomain(p::BlowupMorphism) = p.codomain
center(p::BlowupMorphism) = p.center

function projection(p::BlowupMorphism)
  if !isdefined(p, :projection)
    # create the projection
  end
  return p.projection
end

# TODO: Find better name!
covered_projective_scheme(p::BlowupMorphism) = p.projective_bundle

function exceptional_divisor(p::BlowupMorphism)
  if !isdefined(p, :exceptional_divisor)
    #TODO: Build it!
  end
  return p.exceptional_divisor
end
