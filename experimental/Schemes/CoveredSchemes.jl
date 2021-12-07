export Glueing
export glueing_morphism, patches, glueing_domains, inverse_glueing_morphism

export compose, maximal_extension

mutable struct Glueing{BRT, BRET, RT, RET}
  U::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}	
  V::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}	
  W1::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}	
  W2::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}	
  f::SpecMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}}
  g::SpecMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}}

  function Glueing(U::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}},
      V::Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}},
      f::SpecMor{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}, MPolyPowersOfElement{BRT, BRET, RT, RET}}
    ) where {BRT, BRET, RT, RET}
    is_open_embedding(domain(f), U) || error("the domain of the glueing morphism is not an open subset of the first argument")
    is_open_embedding(codomain(f), V) || error("the codomain of the glueing morphism is not an open subset of the second argument")
    is_isomorphism(f) || error("the given morphism for the glueing is not an isomorphism")
    return new{BRT, BRET, RT, RET}(U, V, domain(f), codomain(f), f, inverse(f))
  end
end

patches(G::Glueing{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET} = [G.U, G.V]
glueing_morphism(G::Glueing{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET} = G.f
inverse_glueing_morphism(G::Glueing{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET} = G.g
glueing_domains(G::Glueing) = [domain(glueing_morphism(G)), codomain(glueing_morphism(G))]
inverse(G::Glueing) = Glueing(G.V, G.U, G.g)

function Base.show(io::IO, G::Glueing)
  print(io, "Glueing of $(patches(G)[1]) and $(patches(G)[2]) along the map $(glueing_morphism(G))")
end

@Markdown.doc """
compose(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}

Given glueings `X ↩ U ↪  Y` and `Y ↩ V ↪ Z`, return the glueing
`X ↩ U ∩ V ↪ Z`. 

**WARNING:** In general such a glueing will not provide a separated scheme. 
Use `maximal_extension` to extend the glueing.
"""
function compose(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}
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
  g = glueing_morphism(H)
  U = preimage(f, domain(g))
  W = preimage(inverse(g), codomain(f))
  V = intersect(codomain(f), domain(g))
  return Glueing(X, Z, compose(restrict(f, U, V), restrict(g, V, W)))
end

@Markdown.doc """
maximal_extension(G::Glueing)

Given a glueing `X ↩ U ↪ Y`, try to find the maximal extension to an open 
subset `V ⊃ U` in both `X` and `Y` so that the resulting scheme is separated.
"""
function maximal_extension(G::Glueing)
  X = patches(G)[1]
  Y = patches(G)[2]
  f = glueing_morphism(G)
  g = inverse_glueing_morphism(G)
  pb_f = pullback(f)
  pb_g = pullback(g)
  U = domain(f)
  V = codomain(f)
  W = localized_ring(OO(X))
  I = ideal(W, one(W))
  for a in images(pb_f)
    D = quotient(ideal(W, denominator(a)) + localized_modulus(OO(X)), ideal(W, numerator(a)))
    I = intersection(I, D)
  end
  extensions = Vector{typeof(G)}()
  for d in gens(saturated_ideal(I))
    e = lifted_denominator(pb_g(OO(U)(d)))
    T = hypersurface_complement(X, d)
    T in [glueing_domains(H)[1] for H in extensions] || (push!(extensions, Glueing(X, Y, SpecMor(T, hypersurface_complement(Y, e), images(pb_f)))))
  end
  return extensions
end

function ==(G::GlueingType, H::GlueingType) where {GlueingType<:Glueing}
  if patches(G)[1] != patches(H)[1]
    return G == inverse(H)
  end
  result = true
  result = result && (patches(G)[2] == patches(H)[2])
  T = glueing_domains(G)[1]
  result = result && (OO(T).(images(pullback(glueing_morphism(G)))) == OO(T).(images(pullback(glueing_morphism(H)))))
  return result
end
