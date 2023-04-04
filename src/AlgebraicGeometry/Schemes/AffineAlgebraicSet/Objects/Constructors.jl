########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    affine_algebraic_set(X::Spec; check::Bool=true) -> AffineAlgebraicSet

Convert `X` to an `AffineAlgebraicSet` by taking the underlying reduced scheme.

If `check=false`, assumes that `X` is already reduced.
"""
function affine_algebraic_set(X::Spec; check::Bool=true)
  if check
    Xred,_ = reduced_scheme(X)
  else
    Xred = X
  end
  AffineAlgebraicSet(Xred, check=check)
end

@doc raw"""
    vanishing_locus(I::MPolyIdeal; check::Bool=true)

Return the vanishing locus of ``I`` as an algebraic set.

This computes the radical of ``I`` if `check=true`
otherwise take on faith that ``I`` is radical.
"""
function vanishing_locus(I::MPolyIdeal{<:MPolyElem}; check::Bool=true)
  if check
    Irad = radical(I)
  else
    Irad = I
  end
  X = Spec(base_ring(Irad), Irad)
  return AffineAlgebraicSet(X, check=check)
end

vanishing_locus(p::MPolyRingElem, check::Bool=true) = vanishing_locus(ideal(parent(p),p), check=check)
########################################################
# (2) Intersections of algebraic sets
########################################################

@doc raw"""
    set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet) -> AbsAffineAlgebraicSet

Return the set theoretic intersection of `X` and `Y` as an AlgebraicSet.

This is the reduced subscheme of the scheme theoretic intersection.
"""
function set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)
  Z = intersect(underlying_scheme(X), underlying_scheme(Y))
  Zred,_ = reduced_scheme(Z)
  # not sure how reduced vs geometrically reduced behaves hence check=true
  return AffineAlgebraicSet(Zred, check=true)
end

########################################################
# (3) Closure of algebraic sets
########################################################

@doc raw"""
  closure(X::AbsAffineAlgebraicSet)

Return the closure of ``X`` in its ambient affine space.
"""
function closure(X::AbsAffineAlgebraicSet)
  Xcl = closure(X, ambient_space(X))
  return affine_algebraic_set(Xcl, check=false)
end
