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
@Markdown.doc """
    compose(G::AbsGlueing, H::AbsGlueing)

Given glueings `X ↩ U ≅ V ↪  Y` and `Y ↩ V' ≅ W ↪ Z`, return the glueing
`X ↩  V ∩ V' ↪ Z`. 

**WARNING:** In general such a glueing will not provide a separated scheme. 
Use `maximal_extension` to extend the glueing.
"""
function compose(G::AbsGlueing, H::AbsGlueing) 
  error("method `compose` not implemented for glueings of type $(typeof(G)) and $(typeof(H))")
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

@Markdown.doc """
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

