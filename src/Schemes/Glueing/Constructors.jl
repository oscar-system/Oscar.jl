export compose, maximal_extension, restrict

########################################################################
# Conversion from SimpleGlueings to Glueings                           #
########################################################################
function Glueing(G::SimpleGlueing)
  X, Y = patches(G)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  incY = inclusion_morphism(V, Y)
  incX = inclusion_morphism(U, X)
  Uo = SpecOpen(U)
  Vo = SpecOpen(V)
  fo = SpecOpenMor(Uo, Vo, [compose(f, incY)], check=false)
  go = SpecOpenMor(Vo, Uo, [compose(g, incX)], check=false)
  return Glueing(X, Y, fo, go, check=false)
end

########################################################################
# Simplified constructor for common special case                       #
########################################################################
function Glueing(
    X::AbsSpec, Y::AbsSpec, 
    f::AbsSpecMor{<:PrincipalOpenSubset}, 
    g::AbsSpecMor{<:PrincipalOpenSubset};
    check::Bool=true)
  return SimpleGlueing(X, Y, f, g, check=check)
end

########################################################################
# Restrictions of Glueings to closed subschemes                        #
########################################################################
function restrict(G::Glueing, X::AbsSpec, Y::AbsSpec; check::Bool=true)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient(U)), ambient(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient(V)), ambient(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  Ures = intersect(X, U)
  Vres = intersect(Y, V)
  return Glueing(X, Y, restrict(f, Ures, Vres, check=check), restrict(g, Vres, Ures, check=check), check=check)
end

function restrict(G::SimpleGlueing, X::AbsSpec, Y::AbsSpec; check::Bool=true)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  UX = PrincipalOpenSubset(X, OO(X)(lifted_numerator(complement_equation(U))))
  VY = PrincipalOpenSubset(Y, OO(Y)(lifted_numerator(complement_equation(V))))
  f_res = restrict(f, UX, VY, check=check)
  g_res = restrict(g, VY, UX, check=check)
  return SimpleGlueing(X, Y, f_res, g_res)
end

########################################################################
# Identification of glueings under isomorphisms of patches             #
########################################################################
@Markdown.doc """
    restrict(G::AbsGlueing, f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)

Given a glueing ``X ↩ U ≅ V ↪ Y`` and isomorphisms ``f : X → X'`` and 
``g: Y → Y'``, return the induced glueing of ``X'`` and ``Y'``.
"""
function restrict(G::AbsGlueing, f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)
  (X1, Y1) = patches(G)
  X1 === domain(f) || error("maps not compatible")
  X2 = codomain(f)
  finv = inverse(f)

  Y1 === domain(g) || error("maps not compatible")
  Y2 = codomain(g)
  ginv = inverse(g)

  (h1, h2) = glueing_morphisms(G)

  U2 = preimage(finv, domain(h1), check=check)
  V2 = preimage(ginv, domain(h2), check=check)

  return Glueing(X2, Y2, 
                 compose(restrict(finv, U2, domain(h1), check=check), 
                         compose(h1, restrict(g, domain(h2), V2, check=check), check=check),
                         check=check
                        ),
                 compose(restrict(ginv, V2, domain(h2), check=check), 
                         compose(h2, restrict(f, domain(h1), U2), check=check), 
                         check=check
                        ),
                 check=check
                )
end
