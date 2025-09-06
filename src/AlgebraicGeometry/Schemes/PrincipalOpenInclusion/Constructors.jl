function complement(X::AbsAffineScheme, h::Vector)
  all(x -> x isa RingElem && parent(x) === OO(X), h) || return complement(X, OO(X).(h))

  U = PrincipalOpenSubset(X, h)
  f = morphism(U, X, gens(OO(X)), check=false)
  inc = PrincipalOpenEmbedding(f, h)
  inc.inverse_on_image = id_hom(U)
  inc.image = U
  return U, inc
end

function complement(X::AbsAffineScheme, h::Any)
  (h isa RingElem && parent(h) === OO(X)) || return complement(X, OO(X)(h))
  return complement(X, [h])
end

