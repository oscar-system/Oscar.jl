
function ==(f::ProjectiveSchemeMor, g::ProjectiveSchemeMor) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for s in gens(homogeneous_coordinate_ring(codomain(f)))
    iszero(pullback(f)(s) - pullback(g)(s)) || return false
  end
  return true
end

function ==(f::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}},
            g::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}})
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return map_on_affine_cones(f) == map_on_affine_cones(g)
end

@doc raw"""
    covered_scheme_morphism(f::ProjectiveSchemeMor)

Given a morphism of `ProjectiveScheme`s ``f : X → Y``, construct and 
return the same morphism as a `CoveredSchemeMorphism` of the `covered_scheme`s 
of ``X`` and ``Y``, respectively.
"""
@attr function covered_scheme_morphism(f::ProjectiveSchemeMor)
  PX = domain(f)
  PY = codomain(f)
  SX = ambient_coordinate_ring(PX)
  SY = ambient_coordinate_ring(PY)
  pbf = pullback(f) # The pullback on the free polynomial rings, not the quotients

  X = covered_scheme(PX)
  Y = covered_scheme(PY)

  mor_dict = IdDict{AbsSpec, AbsSpecMor}()
  U = affine_charts(X)
  for i in 1:ngens(SX)
    U_i = U[i]
    dehom = dehomogenization_map(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
    for j in 1:ngens(SY)
      y = gens(SY, j)
      denom = dehom(pbf(y))
      V_j = affine_charts(Y)[j]
      U_ij = PrincipalOpenSubset(U_i, denom)
      u = inv(OO(U_ij)(denom))
      mor_dict[U_ij] = SpecMor(U_ij, V_j, 
                               hom(OO(V_j), OO(U_ij), 
                                   [OO(U_ij)(dehom(pbf(gens(SY, k))))*u for k in 1:ngens(SY) if k != j]
                                  )
                              )
    end
  end
  # We skip the glueings for the time being.
  # Eventually, they should be made lazy.
  CC = Covering(collect(keys(mor_dict)), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}())
  inherit_glueings!(CC, default_covering(X))
  phi = CoveringMorphism(CC, default_covering(Y), mor_dict)
  push!(coverings(X), CC)

  ff = CoveredSchemeMorphism(X, Y, phi)
  return ff
end

