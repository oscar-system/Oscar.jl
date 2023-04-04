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
    # additional checks in the constructor
    # note that Irad knows that it is a radical ideal
    # by construction
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

########################################################
# (3) Irreducible Components
########################################################

@doc raw"""
    irreducible_components(X::AbsProjectiveAlgebraicSet) -> Vector{ProjectiveVariety}

Return the irreducible components of `X` defined over the same base field.

Note that even if `X` is irreducible, there may be several geometric irreducible components.
"""
function irreducible_components(X::AbsProjectiveAlgebraicSet)
  I = defining_ideal(X)
  J = minimal_primes(I)
  return typeof(X)[vanishing_locus(j,check=false) for j in J]
end


@doc raw"""
    geometric_irreducible_components(X::AbsProjectiveAlgebraicSet) -> Vector{ProjectiveVariety}

Return the geometric irreducible components of `X`.

They are the irreducible components of `X` seen over an algebraically closed field.

This is expensive and involves taking field extensions.
"""
function geometric_irreducible_components(X::AbsProjectiveAlgebraicSet)
  throw(NotImplementedError(:geometric_irreducible_components, X))
end
