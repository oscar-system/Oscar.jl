########################################################################
# Printing                                                             #
########################################################################
function Base.show(io::IO, G::AbsGlueing)
  print(io, "Glueing of $(patches(G)[1]) and $(patches(G)[2]) along the map $(glueing_morphisms(G)[1])")
end

function Base.show(io::IO, G::Glueing)
  print(io, "Glueing of $(patches(G)[1]) and $(patches(G)[2]) along the map $(glueing_morphisms(G)[1])")
end

########################################################################
# Composition of Glueings                                              #
########################################################################
@doc raw"""
    compose(G::AbsGlueing, H::AbsGlueing)

Given glueings `X ↩ U ≅ V ↪  Y` and `Y ↩ V' ≅ W ↪ Z`, return the glueing
`X ↩  V ∩ V' ↪ Z`. 

**WARNING:** In general such a glueing will not provide a separated scheme. 
Use `maximal_extension` to extend the glueing.
"""
function compose(G::AbsGlueing, H::AbsGlueing) 
  # In the default case, we rely on some conversion to the most complicated 
  # glueing and rely on the implementation for composition there.
  return compose(Glueing(G), Glueing(H))
end

function compose(G::Glueing, H::SimpleGlueing) 
  return compose(G, Glueing(H))
end

function compose(G::SimpleGlueing, H::Glueing) 
  return compose(Glueing(G), H)
end

function compose(G::Glueing, H::Glueing) 
  # make sure that Y is the second patch of the first glueing and 
  # the first patch of the second
  if patches(G)[2] === patches(H)[2]
    return compose(G, inverse(H))
  elseif patches(G)[1] === patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] === patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y === patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = preimage(f, domain(g))
  W_new = preimage(g_inv, codomain(f))
  V_new = intersect(codomain(f), domain(g))
  return Glueing(X, Z, 
             compose(restrict(f, U_new, V_new), restrict(g, V_new, W_new)),
             compose(restrict(g_inv, W_new, V_new), restrict(f_inv, V_new, U_new))
	     )
end

@doc raw"""
    maximal_extension(G::Glueing)

Given a glueing `X ↩ U ≅ V ↪ Y`, try to find the maximal extension to an open 
subset `U' ⊃ U` in `X` and `V' ⊃ V` in `Y` so that the resulting scheme is separated.
"""
function maximal_extension(G::Glueing)
  X = patches(G)[1]
  Y = patches(G)[2]
  f, g = glueing_morphisms(G)
  f_ext = maximal_extension(X, Y, generic_fractions(f))
  g_ext = maximal_extension(Y, X, generic_fractions(g))
  U_new = preimage(f_ext, domain(g_ext))
  V_new = preimage(g_ext, U_new)
  issubset(domain(g_ext), V_new) || error("extension failed")
  f_ext = restrict(f_ext, U_new, V_new)
  g_ext = restrict(g_ext, V_new, U_new)
  return Glueing(X, Y, f_ext, g_ext)
end

function compose(G::SimpleGlueing, H::SimpleGlueing)
  if patches(G)[2] === patches(H)[2]
    return compose(G, inverse(H))
  elseif patches(G)[1] === patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] === patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y === patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = PrincipalOpenSubset(ambient_scheme(domain(f)), 
                              [complement_equation(domain(f)), 
                               lifted_numerator(pullback(f)(complement_equation(domain(g))))
                              ])
  W_new = PrincipalOpenSubset(ambient_scheme(domain(g_inv)), 
                              [complement_equation(domain(g_inv)), 
                               lifted_numerator(pullback(g_inv)(complement_equation(domain(f_inv))))
                              ])
  V_new = PrincipalOpenSubset(ambient_scheme(domain(g)), 
                              [complement_equation(domain(g)), complement_equation(domain(f_inv))]
                             )
  h = compose(restrict(f, U_new, V_new, check=false), 
              restrict(g, V_new, W_new, check=false))
  h_inv = compose(restrict(g_inv, W_new, V_new, check=false),
                  restrict(f_inv, V_new, U_new, check=false))
  set_attribute!(h, :inverse, h_inv)
  set_attribute!(h_inv, :inverse, h)
  return SimpleGlueing(X, Z, h, h_inv)
end

########################################################################
# Equality tests                                                       #
########################################################################
function ==(G::AbsGlueing, H::AbsGlueing)
  if patches(G)[1] != patches(H)[1]
    return G == inverse(H)
  end
  patches(G)[2] === patches(H)[2] || return false
  glueing_morphisms(G) == glueing_morphisms(H) || return false
  return true
end

########################################################################
# Base change
########################################################################
struct BaseChangeGlueingData{T1, T2}
  phi::T1
  G::T2
  patch_change1::AbsSpecMor
  patch_change2::AbsSpecMor
end

function _compute_glueing_base_change(gd::BaseChangeGlueingData)
  # extraction of the glueing data
  G = gd.G
  pc1 = gd.patch_change1
  pc2 = gd.patch_change2
  phi = gd.phi

  # computation of the new glueing
  (X, Y) = patches(G)
  (U, V) = glueing_domains(G)
  (f, g) = glueing_morphisms(G)
  
  XX = domain(pc1)
  YY = domain(pc2)

  UU, map_U = base_change(phi, U, ambient_map=pc1)
  VV, map_V = base_change(phi, V, ambient_map=pc2)

  # TODO: The following will not yet work for glueings along SpecOpens.
  U isa PrincipalOpenSubset && V isa PrincipalOpenSubset || error("base change not implemented for glueings along SpecOpens")

  _, ff, _ = base_change(phi, f, domain_map=map_U, codomain_map=map_V)
  _, gg, _ = base_change(phi, g, domain_map=map_V, codomain_map=map_U)

  return Glueing(XX, YY, ff, gg)
end

function base_change(phi::Any, G::AbsGlueing;
    patch_change1::AbsSpecMor=base_change(phi, glueing_domains(G)[1])[2], # The base change morphism for 
    patch_change2::AbsSpecMor=base_change(phi, glueing_domains(G)[2])[2]  # the two glueing patches
  )
  gd = BaseChangeGlueingData{typeof(phi), typeof(G)}(phi, G, patch_change1, patch_change2)
  return LazyGlueing(domain(patch_change1), domain(patch_change2), _compute_glueing_base_change, gd)
end

