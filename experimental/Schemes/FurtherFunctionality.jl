@attr Bool function is_smooth(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  R = base_ring(OO(X))
  L = localized_ring(OO(X))
  I = localized_modulus(OO(X))
  f = gens(I)
  Df = jacobi_matrix(f)

  A = map_entries(x->OO(X)(x), Df)
  success, _ = Oscar._is_projective(A, X)
  return success
end

is_smooth(X::AbsSpec{<:Ring, <:MPolyRing}) = true
is_smooth(X::AbsSpec{<:Ring, <:MPolyLocalizedRing}) = true

@attr Bool function is_smooth(X::AbsSpec{<:Ring, <:MPolyQuo})
  R = base_ring(OO(X))
  I = modulus(OO(X))
  f = gens(I)
  Df = jacobi_matrix(f)

  A = map_entries(x->OO(X)(x), Df)
  success, _ = Oscar._is_projective(A, X)
  return success
end

