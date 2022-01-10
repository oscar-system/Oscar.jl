export Glueing
export glueing_morphism, patches, glueing_domains, inverse_glueing_morphism, inverse

export compose, maximal_extension

export CoveredScheme

import Oscar.Graphs: Graph, Directed

@Markdown.doc """
Glueing{BRT, BRET, RT, RET}

Glueing of two affine schemes ``X ↩ U ≅ V ↪ Y`` along open subsets 
``U ⊂ X`` and ``V ⊂ Y via some isomorphism ``φ : U → V``.
"""
mutable struct Glueing{BRT, BRET, RT, RET}
  X::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}	
  Y::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}	
  U::SpecOpen{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}
  V::SpecOpen{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}
  f::SpecOpenMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}}
  g::SpecOpenMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}}

  function Glueing(X::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}},
      Y::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}},
      f::SpecOpenMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}},
      g::SpecOpenMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}}
    ) where {BRT, BRET, RT, RET}
    parent(domain(f)) == X || error("the domain of the glueing morphism is not an open subset of the first argument")
    parent(codomain(f)) == Y || error("the codomain of the glueing morphism is not an open subset of the second argument")
    (domain(f) == codomain(g) && domain(g) == codomain(f)) || error("maps can not be isomorphisms")
    return new{BRT, BRET, RT, RET}(X, Y, domain(f), codomain(f), f, g)
  end
end

patches(G::Glueing{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET} = G.X, G.Y
glueing_morphism(G::Glueing{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET} = G.f
inverse_glueing_morphism(G::Glueing{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET} = G.g
glueing_domains(G::Glueing) = domain(glueing_morphism(G)), codomain(glueing_morphism(G))
inverse(G::Glueing) = Glueing(G.Y, G.X, G.g, G.f)

function Base.show(io::IO, G::Glueing)
  print(io, "Glueing of $(patches(G)[1]) and $(patches(G)[2]) along the map $(glueing_morphism(G))")
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
  X = patches(G)[1]
  Y = patches(G)[2]
  Y == patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f = glueing_morphism(G)
  f_inv = inverse_glueing_morphism(G)
  g = glueing_morphism(H)
  g_inv = inverse_glueing_morphism(H)
  U_new = preimage(f, domain(g))
  W_new = preimage(g_inv, codomain(f))
  V_new = intersect(codomain(f), domain(g))
  return Glueing(X, Z, 
             compose(restriction(f, U_new, V_new), restriction(g, V_new, W_new)),
             compose(restriction(g_inv, W_new, V_new), restriction(f_inv, V_new, U_new))
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
  f = glueing_morphism(G)
  g = inverse_glueing_morphism(G)
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
  glueing_morphism(G) == glueing_morphism(H) || return false
  return true
end


mutable struct CoveredScheme{BRT, BRET, RT, RET} <:Scheme{BRT, BRET}
  patches::Vector{Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}}
  glueings::Dict{Tuple{Spec{BRT, BRET, RT, RET}, Spec{BRT, BRET, RT, RET}}, Glueing{BRT, BRET, RT, RET}}

  # fields for caching
  glueing_graph::Graph{Directed}

  function CoveredScheme(
      patches::Vector{Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}},
      glueings::Dict{Tuple{Spec{BRT, BRET, RT, RET}, Spec{BRT, BRET, RT, RET}}, Glueing{BRT, BRET, RT, RET}}
    ) where {BRT, BRET, RT, RET}
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(base_ring(OO(patches[1])))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(g, (Y, X))
	inverse(g[(X, Y)]) == g[(Y, X)] || error("glueings are not inverse of each other")
      else
	g[(Y, X)] = inverse(g[(X, Y)])
      end
    end
    return new{BRT, BRET, RT, RET}(X, glueings)
  end
end

patches(X::CoveredScheme) = X.patches
npatches(X::CoveredScheme) = length(X.patches)
glueings(X::CoveredScheme) = X.glueings
getindex(X::CoveredScheme, i::Int) = X.patches[i]

function glueing_graph(X::CoveredScheme) 
  if !isdefined(X, :glueing_graph)
    n = npatches(X)
    gg = Graph{Directed}(n)
    for (i,j) in keys(glueings(X))
      g = X[(i,j)] 
      Y = X[i]
      Z = X[j]
      (U, V) = glueing_domains(g)
      is_dense(U, Y) && add_edge!(gg, i, j)
      is_dense(V, Z) && add_edge!(gg, j, i)
    end
    X.glueing_graph = gg
  end
  return X.glueing_graph
end

function are_glued(X::CoveredScheme, i::Int64, j::Int64)
  return exists_path(glueing_graph(X), i, j)
end

function getindex(X::CoveredScheme, a::Tuple{Int64, Int64})
  if haskey(glueings(X), a)
    return glueings(X)[a]
  end
  if are_glued(X, a[1], a[2])
    p = get_path(glueing_graph(X), a[1], a[2])
    # construct the glueing along p
    glueings(X)[a] = g
    return g
  end
  glueings(X)[a] = trivial_glueing(X[a[1]], X[a[2]])
  glueings(X)[(a[2], a[1])] = inverse(glueings(X)[a])
  return glueings(X)[a]
end
    
