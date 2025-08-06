########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    algebraic_set(X::AbsProjectiveScheme; is_reduced::Bool=false, check::Bool=true) -> ProjectiveAlgebraicSet

Convert `X` to a `ProjectiveAlgebraicSet` by considering its underlying reduced scheme.

If `is_reduced` is `true` assume that `X` is already reduced.

```jldoctest
julia> P, (x0, x1, x2) = graded_polynomial_ring(QQ,[:x0,:x1,:x2]);

julia> X = proj(ideal([x0*x1^2, x2]))
Projective scheme
  over rational field
defined by ideal (x0*x1^2, x2)

julia> Y = algebraic_set(X)
Projective algebraic set
  in projective 2-space over QQ with coordinates [x0, x1, x2]
defined by ideal (x0*x1^2, x2)

```
"""
function algebraic_set(X::AbsProjectiveScheme; is_reduced::Bool=false, check::Bool=true)
  return ProjectiveAlgebraicSet(X, is_reduced=is_reduced, check=check)
end

@doc raw"""
    algebraic_set(I::MPolyIdeal{MPolyDecRingElem})

Return the projrective algebraic set defined by the homogeneous ideal ``I``.

```jldoctest
julia> P,(x0,x1) = graded_polynomial_ring(QQ,[:x0,:x1]);

julia> algebraic_set(ideal([x0,x1]))
Projective algebraic set
  in projective 1-space over QQ with coordinates [x0, x1]
defined by ideal (x0, x1)

```
"""
function algebraic_set(I::MPolyIdeal{<:MPolyDecRingElem};
                       is_radical::Bool=false,
                       check::Bool=true)
  X = proj(base_ring(I), I)
  return algebraic_set(X, is_reduced=is_radical, check=check)
end


@doc raw"""
    algebraic_set(p::MPolyDecRingElem; check::Bool=true)

Return the projective algebraic set defined by the homogeneous polynomial `p`.
"""
algebraic_set(p::MPolyDecRingElem; is_radical::Bool=false, check::Bool=true) = algebraic_set(ideal(parent(p),p), is_radical=is_radical, check=check)

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
  Z = intersect(fat_scheme(X), fat_scheme(Y))
  # not sure how reduced vs geometrically reduced behaves hence check=true
  return algebraic_set(Z)
end

########################################################
# (3) Irreducible Components
########################################################

@doc raw"""
    irreducible_components(X::AbsProjectiveAlgebraicSet) -> Vector{ProjectiveVariety}

Return the irreducible components of ``X`` defined over the base field of ``X``.

Note that even if ``X`` is irreducible, there may be several geometrically irreducible components.

```jldoctest
julia> P1 = projective_space(QQ,1)
Projective space of dimension 1
  over rational field
with homogeneous coordinates [s0, s1]

julia> (s0,s1) = homogeneous_coordinates(P1);

julia> X = algebraic_set((s0^2+s1^2)*s1)
Projective algebraic set
  in projective 1-space over QQ with coordinates [s0, s1]
defined by ideal (s0^2*s1 + s1^3)

julia> (X1,X2) = irreducible_components(X)
2-element Vector{ProjectiveAlgebraicSet{QQField, MPolyQuoRing{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 V(s0^2 + s1^2)
 V(s1)

julia> X1  # irreducible but not geometrically irreducible
Projective algebraic set
  in projective 1-space over QQ with coordinates [s0, s1]
defined by ideal (s0^2 + s1^2)

```
"""
function irreducible_components(X::AbsProjectiveAlgebraicSet{S,T}) where {S,T}
  I = fat_ideal(X)
  J = minimal_primes(I)
  return ProjectiveAlgebraicSet{S,T}[algebraic_set(j, is_radical=true, check=false) for j in J]
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

Base.union(X::AbsProjectiveAlgebraicSet, Y::AbsProjectiveAlgebraicSet) = ProjectiveAlgebraicSet(union(fat_scheme(X), fat_scheme(Y)))

