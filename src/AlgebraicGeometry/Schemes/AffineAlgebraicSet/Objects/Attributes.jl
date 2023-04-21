########################################################################
#
# (1) AbsSpec interface
#
########################################################################

@doc raw"""
    underlying_scheme(X::AffineAlgebraicSet) -> AbsSpec

Return the underlying reduced scheme defining ``X``.

This is used to forward the `AbsSpec` functionality to ``X``, but may
trigger the computation of a radical ideal. Hence this can be expensive.
"""
function underlying_scheme(X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  if isdefined(X, :Xred)
    return X.Xred
  end
  # avoid constructing a morphism
  I = ambient_closure_ideal(overlying_scheme(X))
  Irad = radical(I)
  X.Xred = Spec(base_ring(Irad), Irad)
  return X.Xred
end

function underlying_scheme(X::AffineAlgebraicSet)
  if isdefined(X, :Xred)
    return X.Xred
  end
  X.Xred = reduced_scheme(X.X)[1]
  return X.Xred
end
########################################################################
#
# (2) AbsAffineAlgebraicSet interface
#
########################################################################
@doc raw"""
    overlying_scheme(X::AffineAlgebraicSet) -> AbsSpec

Return a scheme whose reduced subscheme is ``X``.

This does not trigger any computation and is therefore cheap.
Use this instead of `underlying_scheme` when possible.
"""
function overlying_scheme(X::AffineAlgebraicSet)
  # if we know the reduced structure already, we can return that.
  if isdefined(X, :Xred)
    return X.Xred
  end
  return X.X
end


########################################################################
#
# (3) Further attributes
#
########################################################################

@doc raw"""
    vanishing_ideal(X::AbsAffineAlgebraicSet) -> Ideal

Return the radical ideal of all polynomials vanishing in ``X``.

!!! note
    This involves the computation of a radical which is expensive.
"""
vanishing_ideal(X::AbsAffineAlgebraicSet) = ambient_closure_ideal(X)

@doc raw"""
    ideal(X::AbsAffineAlgebraicSet) -> Ideal

Return an ideal whose radical is the vanishing ideal of `X`.
"""
ideal(X::AbsAffineAlgebraicSet) = ambient_closure_ideal(overlying_scheme(X))


# avoid computing the underlying scheme
ambient_space(X::AbsAffineAlgebraicSet) = ambient_space(overlying_scheme(X))



