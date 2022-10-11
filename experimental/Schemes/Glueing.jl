export glueing_morphisms, patches, glueing_domains, inverse_glueing_morphism, inverse
export glueing_type

export compose, maximal_extension, restriction

########################################################################
# The interface for abstract glueings of affine schemes                #
########################################################################
function patches(G::AbsGlueing)
  return patches(underlying_glueing(G))
end

function glueing_morphisms(G::AbsGlueing)
  return glueing_morphisms(underlying_glueing(G))
end

function glueing_domains(G::AbsGlueing)
  return glueing_domains(underlying_glueing(G))
end

@attr function inverse(G::AbsGlueing)
  Ginv = G(inverse(underlying_glueing(G)))
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end

function Base.show(io::IO, G::AbsGlueing)
  print(io, "Glueing of $(patches(G)[1]) and $(patches(G)[2]) along the map $(glueing_morphisms(G)[1])")
end

########################################################################
# Methods for Glueing                                                  #
########################################################################
patches(G::Glueing) = G.X, G.Y
glueing_morphisms(G::Glueing) = G.f, G.g
glueing_domains(G::Glueing) = domain(G.f), domain(G.g)
inverse(G::Glueing) = Glueing(G.Y, G.X, G.g, G.f, check=false)

function Base.show(io::IO, G::Glueing)
  print(io, "Glueing of $(patches(G)[1]) and $(patches(G)[2]) along the map $(glueing_morphisms(G)[1])")
end

@Markdown.doc """
compose(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}

Given glueings `X ↩ U ≅ V ↪  Y` and `Y ↩ V' ≅ W ↪ Z`, return the glueing
`X ↩  V ∩ V' ↪ Z`. 

**WARNING:** In general such a glueing will not provide a separated scheme. 
Use `maximal_extension` to extend the glueing.
"""
function compose(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}
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

# A simplified constructor for the common case to glue along principal open sets
function Glueing(
    X::SpecType, p::RET, Y::SpecType, q::RET,
    a::Vector{FracType}, b::Vector{FracType};
    check::Bool=true
  ) where {
           SpecType<:Spec, 
           RET<:MPolyElem,
           FracType
          }
  parent(p) == base_ring(OO(X)) || error("polynomial does not belong to the correct ring")
  parent(q) == base_ring(OO(Y)) || error("polynomial does not belong to the correct ring")
  UU = SpecOpen(X, [p])
  VV = SpecOpen(Y, [q])
  U = UU[1]
  V = VV[1]
  f = SpecMor(U, Y, a)
  g = SpecMor(V, X, b)
  if check
    compose(restrict(f,U, V), restrict(g, V, U)) == identity_map(U) || error("glueing maps are not inverse to each other")
    compose(restrict(g, V, U), restrict(f, U, V)) == identity_map(V) || error("glueing maps are not inverse to each other")
  end
  return Glueing(X, Y, SpecOpenMor(UU, VV, [f]), SpecOpenMor(VV, UU, [g]))
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

function ==(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}
  if patches(G)[1] != patches(H)[1]
    return G == inverse(H)
  end
  patches(G)[2] === patches(H)[2] || return false
  glueing_morphisms(G) == glueing_morphisms(H) || return false
  return true
end

function glueing_type(X::SpecType) where {SpecType<:Spec}
  return Glueing{SpecType, open_subset_type(SpecType), morphism_type(open_subset_type(SpecType), open_subset_type(SpecType))}
end

function glueing_type(::Type{SpecType}) where {SpecType<:Spec}
  return Glueing{SpecType, open_subset_type(SpecType), morphism_type(open_subset_type(SpecType), open_subset_type(SpecType))}
end

function restrict(G::Glueing, X::SpecType, Y::SpecType; check::Bool=true) where {SpecType<:Spec}
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient(U)), ambient(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient(V)), ambient(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  return Glueing(X, Y, restrict(f, X, Y, check=check), restrict(g, Y, X, check=check), check=check)
end

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
