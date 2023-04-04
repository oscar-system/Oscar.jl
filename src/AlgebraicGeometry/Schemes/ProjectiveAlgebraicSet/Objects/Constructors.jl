########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true) -> ProjectiveAlgebraicSet

Convert `X` to an `ProjectiveAlgebraicSet` by taking the underlying reduced scheme.

If `check=false`, assumes that `X` is already reduced.
"""
function projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true)
  if check
    Xred,_ = reduced_scheme(X)
  else
    Xred = X
  end
  ProjectiveAlgebraicSet(Xred, check=check)
end

@doc raw"""
    vanishing_locus(I::MPolyIdeal{MPolyDecRingElem}; check::Bool=true)

Return the vanishing locus of ``I`` as an algebraic set in projective space.

This computes the radical of ``I`` if `check=true`
otherwise take on faith that ``I`` is radical.
"""
function vanishing_locus(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
  if check
    Irad = radical(I)
  else
    Irad = I
  end
  X = ProjectiveScheme(base_ring(Irad),Irad)
  return ProjectiveAlgebraicSet(X, check=check)
end

vanishing_locus(p::MPolyDecRingElem; check::Bool=true) = vanishing_locus(ideal(parent(p),p), check=check)
########################################################
# (2) Intersections of algebraic sets
########################################################

@doc raw"""
    set_theoretic_intersection(X::AbsProjectiveAlgebraicSet, Y::AbsProjectiveAlgebraicSet) -> AbsProjectiveAlgebraicSet

Return the set theoretic intersection of `X` and `Y` as a `ProjectiveAlgebraicSet`.

This is the reduced subscheme of the scheme theoretic intersection.
"""
function set_theoretic_intersection(X::AbsProjectiveAlgebraicSet, Y::AbsProjectiveAlgebraicSet)
  throw(NotImplementedError(:set_theoretic_intersection,X))
  Z = intersect(underlying_scheme(X), underlying_scheme(Y))
  Zred,_ = reduced_scheme(Z)
  # not sure how reduced vs geometrically reduced behaves hence check=true
  return ProjectiveAlgebraicSet(Zred, check=true)
end
