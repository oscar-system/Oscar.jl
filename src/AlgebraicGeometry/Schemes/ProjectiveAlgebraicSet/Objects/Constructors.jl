########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true) -> ProjectiveAlgebraicSet

Convert `X` to an `ProjectiveAlgebraicSet` by taking the underlying reduced scheme.

If `check=false`, assumes that `X` is already reduced.

```jldoctest
julia> P,(x0,x1,x2) = graded_polynomial_ring(QQ,[:x0,:x1,:x2]);

julia> X = projective_scheme(ideal([x0*x1^2, x2]))
Projective scheme
  over Rational Field
  defined by
ideal(x0*x1^2, x2)

julia> Y = projective_algebraic_set(X)
Vanishing locus
  in IP^2 over Rational Field
  of ideal(x2, x0*x1)

```
"""
function projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true)
  if check
    Xred = reduced_scheme(X)
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

```jldoctest
julia> P,(x0,x1) = graded_polynomial_ring(QQ,[:x0,:x1]);

julia> vanishing_locus(ideal([x0,x1]))
Vanishing locus
  in IP^1 over Rational Field
  of ideal(x1, x0)

```
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

```jldoctest
julia> P1 = projective_space(QQ,1)
Projective space of dimension 1
  over Rational Field

julia> (s0,s1) = homogeneous_coordinates(P1);

julia> X = vanishing_locus((s0^2+s1^2)*s1)
Vanishing locus
  in IP^1 over Rational Field
  of ideal(s0^2*s1 + s1^3)

julia> (X1,X2) = irreducible_components(X)
2-element Vector{ProjectiveAlgebraicSet{QQField, MPolyQuoRing{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 Projective algebraic set in IP^1 over Rational Field
 Projective algebraic set in IP^1 over Rational Field

julia> X1  # irreducible but not geometrically irreducible
Vanishing locus
  in IP^1 over Rational Field
  of ideal(s0^2 + s1^2)

```
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
