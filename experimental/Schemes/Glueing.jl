export Glueing
export glueing_morphisms, patches, glueing_domains, inverse_glueing_morphism, inverse
export glueing_type

export compose, maximal_extension, restriction


@Markdown.doc """
Glueing{SpecType<:Spec, OpenType<:SpecOpen, MorType<:SpecOpenMor}

Glueing of two affine schemes ``X ↩ U ≅ V ↪ Y`` along open subsets 
``U ⊂ X`` and ``V ⊂ Y via some isomorphism ``φ : U → V``.
"""
mutable struct Glueing{SpecType<:Spec, OpenType<:SpecOpen, MorType<:SpecOpenMor}
  X::SpecType	
  Y::SpecType
  U::OpenType
  V::OpenType
  f::MorType # f : U → V 
  g::MorType

  function Glueing(
      X::SpecType, Y::SpecType, f::MorType, g::MorType; check::Bool=true
    ) where {
      SpecType<:Spec, MorType<:SpecOpenMor
    }
    ambient(domain(f)) === X || error("the domain of the glueing morphism is not an open subset of the first argument")
    ambient(codomain(f)) === Y || error("the codomain of the glueing morphism is not an open subset of the second argument")
    if check
      (is_canonically_isomorphic(domain(f), codomain(g)) && 
       is_canonically_isomorphic(domain(g), codomain(f))) || error("maps can not be isomorphisms")
      compose(f, g) == identity_map(domain(f)) || error("glueing maps are not inverse of each other")
      compose(g, f) == identity_map(domain(g)) || error("glueing maps are not inverse of each other")
    end
    return new{SpecType, open_subset_type(X), MorType}(X, Y, domain(f), domain(g), f, g)
  end
end

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
  if patches(G)[2] == patches(H)[2] 
    return compose(G, inverse(H))
  elseif patches(G)[1] == patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] == patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y == patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = preimage(f, domain(g))
  W_new = preimage(g_inv, codomain(f))
  V_new = intersect(codomain(f), domain(g))
  return Glueing(X, Z, 
             compose(restriction(f, U_new, V_new), restriction(g, V_new, W_new)),
             compose(restriction(g_inv, W_new, V_new), restriction(f_inv, V_new, U_new))
	     )
end

# A simplified constructor for the common case to glue along principal open sets
function Glueing(
    X::SpecType, p::RET, Y::SpecType, q::RET,
    a::Vector{FracType}, b::Vector{FracType};
    check::Bool=true
  ) where {
           SpecType<:Spec, 
           SpecMorType<:SpecMor,
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
  f_ext = restriction(f_ext, preimage(f_ext, domain(g_ext)), domain(g_ext))
  g_ext = restriction(g_ext, preimage(g_ext, domain(f_ext)), domain(f_ext))
  return Glueing(X, Y, f_ext, g_ext)
end

function ==(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}
  if patches(G)[1] != patches(H)[1]
    return G == inverse(H)
  end
  patches(G)[2] == patches(H)[2] || return false
  glueing_morphisms(G) == glueing_morphisms(H) || return false
  return true
end

function glueing_type(X::SpecType) where {SpecType<:Spec}
  return Glueing{SpecType, open_subset_type(SpecType), morphism_type(open_subset_type(SpecType), open_subset_type(SpecType))}
end

function glueing_type(::Type{SpecType}) where {SpecType<:Spec}
  return Glueing{SpecType, open_subset_type(SpecType), morphism_type(open_subset_type(SpecType), open_subset_type(SpecType))}
end

function restriction(G::Glueing, X::SpecType, Y::SpecType; check::Bool=true) where {SpecType<:Spec}
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient(U)), ambient(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient(V)), ambient(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  return Glueing(X, Y, restriction(f, X, Y, check=check), restriction(g, Y, X, check=check), check=check)
end
