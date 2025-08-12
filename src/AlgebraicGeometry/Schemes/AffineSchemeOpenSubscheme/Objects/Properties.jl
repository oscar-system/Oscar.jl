
########################################################################
# Properties of AffineSchemeOpenSubschemes                                              #
########################################################################

function is_dense(U::AffineSchemeOpenSubscheme)
  X = ambient_scheme(U)
  I = [modulus(OO(closure(V, X))) for V in affine_patches(U)]
  J = pre_image_ideal(ideal(OO(X), [one(OO(X))]))
  for i in I
    J = intersect(J, i)
  end
  return J == modulus(OO(X))
end


function is_dense(U::AffineSchemeOpenSubscheme{<:AbsAffineScheme{BRT,<:MPolyQuoRing}}) where {BRT}
  X = ambient_scheme(U)
  I = [modulus(OO(closure(V, X))) for V in affine_patches(U)]
  R = ambient_coordinate_ring(X)
  J = ideal(R,[R(1)])
  for i in I
    J = intersect(J, i)
  end
  return J == modulus(OO(X))
end

function is_dense(U::AffineSchemeOpenSubscheme{<:AbsAffineScheme{BRT, <:Union{MPolyRing,MPolyLocRing}}, BRT}) where {BRT}
  X = ambient_scheme(U)
  return any(V -> closure(V,X)==X, affine_patches(U))
end
