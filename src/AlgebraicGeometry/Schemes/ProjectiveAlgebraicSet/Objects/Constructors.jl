########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true) -> ProjectiveAlgebraicSet

Convert `X` to an `ProjectiveAlgebraicSet` by taking the underlying reduced scheme.

If `check=false`, assumes that `X` is already reduced.

```
julia> P,(x0,x1,x2) = graded_polynomial_ring(QQ,[:x0,:x1,:x2]);

julia> X = projective_scheme(ideal([x0*x1^2, x2]))
Projective scheme
  over Rational field
  defined by ideal(x0*x1^2, x2)

julia> Y = Oscar.projective_algebraic_set(X)
Vanishing locus
  in Projective 2-space over QQ
  of ideal(x2, x0*x1)

```
"""
function projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true)
  @check
  if check
    Xred = reduced_scheme(X)
  else
    Xred = X
  end
  ProjectiveAlgebraicSet(Xred, check=check)
end

@doc raw"""
    vanishing_locus(I::MPolyIdeal{MPolyDecRingElem}; check::Bool=true)

Return the vanishing locus of the homogeneous ideal ``I`` as an algebraic set
in projective space.

This computes the radical of ``I`` if `check=true`.
Otherwise Oscar takes on faith that ``I`` is radical.

```
julia> P,(x0,x1) = graded_polynomial_ring(QQ,[:x0,:x1]);

julia> Oscar.vanishing_locus(ideal([x0,x1]))
Vanishing locus
  in Projective 1-space over QQ
  of ideal(x1, x0)

```
"""
function vanishing_locus(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
  @check
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

projective_algebraic_set(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true) = vanishing_locus(I, check=check)

@doc raw"""
    vanishing_locus(p::MPolyDecRingElem; check::Bool=true)

Return the vanishing locus of the homogeneous polynomial `p` as an algebraic set
in projective space.

This computes the radical of ``I`` if `check=true`
otherwise take on faith that ``I`` is radical.
"""
vanishing_locus(p::MPolyDecRingElem; check::Bool=true) = vanishing_locus(ideal(parent(p),p), check=check)

projective_algebraic_set(p::MPolyDecRingElem; check::Bool=true) = vanishing_locus(p, check=check)
########################################################
# (2) Intersections of algebraic sets
########################################################

@doc raw"""
    set_theoretic_intersection(X::AbsProjectiveAlgebraicSet, Y::AbsProjectiveAlgebraicSet) -> AbsProjectiveAlgebraicSet

Return the set theoretic intersection of `X` and `Y` as as algebraic sets
in projective space.

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

Return the irreducible components of ``X`` defined over the base field of ``X``.

Note that even if ``X`` is irreducible, there may be several geometrically irreducible components.

```
julia> P1 = projective_space(QQ,1)
Projective space of dimension 1
  with homogeneous coordinates s0 s1
  over Rational field

julia> (s0,s1) = homogeneous_coordinates(P1);

julia> X = Oscar.vanishing_locus((s0^2+s1^2)*s1)
Vanishing locus
  in Projective 1-space over QQ
  of ideal(s0^2*s1 + s1^3)

julia> (X1,X2) = Oscar.irreducible_components(X)
2-element Vector{ProjectiveAlgebraicSet{QQField, MPolyQuoRing{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 Vanishing locus in IP^1 of ideal(s0^2 + s1^2)
 Vanishing locus in IP^1 of ideal(s1)

julia> X1  # irreducible but not geometrically irreducible
Vanishing locus
  in Projective 1-space over QQ
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

Return the geometrically irreducible components of `X`.

They are the irreducible components of `X` seen over an algebraically closed field.

This is expensive and involves taking field extensions.
"""
function geometric_irreducible_components(X::AbsProjectiveAlgebraicSet)
  throw(NotImplementedError(:geometric_irreducible_components, X))
end
