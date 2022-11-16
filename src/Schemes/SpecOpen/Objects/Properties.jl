export is_dense

########################################################################
# Properties of SpecOpens                                              #
########################################################################

function is_dense(U::SpecOpen)
  X = ambient_scheme(U)
  I = [modulus(OO(closure(V, X))) for V in affine_patches(U)]
  J = pre_image_ideal(ideal(OO(X), [one(OO(X))]))
  for i in I
    J = intersect(J, i)
  end
  return J == modulus(OO(X))
end

