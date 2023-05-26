########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    affine_algebraic_set(X::Spec; is_reduced=false, check=true) -> AffineAlgebraicSet

Convert the `X` to an `AffineAlgebraicSet` by considering its reduced structure.

If `is_reduced` is set, assume that `X` is already reduced.
If `is_reduced` and `check` are set, check that `X` is actually reduced as claimed.
"""
function affine_algebraic_set(X::Spec; is_reduced::Bool=false, check::Bool=true)
  return AffineAlgebraicSet(X, is_reduced=is_reduced, check=check)
end

@doc raw"""
    vanishing_locus(I::MPolyIdeal; is_radical=false, check::Bool=true)

Return the vanishing locus of ``I`` as an affine algebraic set.

If it is known a priori that `I` is a radical ideal and the base field perfect,
one can set `is_radical=true` and `check=false` to speed up some
computations.

```jldoctest
julia> R, (x,y) = GF(2)[:x,:y];

julia> X = vanishing_locus(ideal([y^2+y+x^3+1,x]))
Vanishing locus
  in Affine 2-space over GF(2)
  of ideal(x, y^2 + y + 1)

```
"""
function vanishing_locus(I::MPolyIdeal{<:MPolyElem}; is_radical=false, check::Bool=true)
  X = Spec(base_ring(I), I)
  return AffineAlgebraicSet(X, is_reduced=is_radical, check=check)
end

affine_algebraic_set(I::MPolyIdeal{<:MPolyElem}; kwargs...) = vanishing_locus(I; kwargs...)

@doc raw"""
    vanishing_locus(p::MPolyRingElem; is_radical=false, check::Bool=true)

Return the vanishing locus of the multivariate polynomial `p`
seen as an affine algebraic set.

If it is known a priori that `p` is radical and the base field perfect, one can set
`is_radical=true` and `check=false` to speed up some computations.

```jldoctest
julia> R, (x,y) = QQ[:x,:y];

julia> X = vanishing_locus((y^2+y+x^3+1)*x^2)
Vanishing locus
  in Affine 2-space over QQ
  of ideal(x^4 + x*y^2 + x*y + x)

julia> R, (x,y) = GF(2)[:x,:y];

julia> X = vanishing_locus((y^2+y+x^3+1)*x^2)
Vanishing locus
  in Affine 2-space over GF(2)
  of ideal(x^4 + x*y^2 + x*y + x)

```
"""
vanishing_locus(p::MPolyRingElem, check::Bool=true) = vanishing_locus(ideal(parent(p),p), check=check)

affine_algebraic_set(p::MPolyElem, check::Bool) = vanishing_locus(p, check=check)

########################################################
# (2) Intersections of algebraic sets
########################################################

@doc raw"""
    set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)

Return the set theoretic intersection of `X` and `Y` as an algebraic set.

```jldoctest set_theoretic_intersection
julia> A = affine_space(QQ, [:x,:y])
Affine space of dimension 2
  with coordinates x y
  over Rational Field

julia> X = vanishing_locus(ideal([y - x^2]))
Vanishing locus
  in Affine 2-space over Rational Field
  of ideal(-x^2 + y)

julia> Y = vanishing_locus(ideal([y]))
Vanishing locus
  in Affine 2-space over Rational Field
  of ideal(y)

julia> Zred = set_theoretic_intersection(X, Y)
Vanishing locus
  in Affine 2-space over Rational Field
  of ideal(-x^2 + y, y)


```
Note that the set theoretic intersection forgets the intersection multiplicities
which the scheme theoretic intersection remembers. Therefore they are different.

```jldoctest set_theoretic_intersection
julia> vanishing_ideal(Zred); # computes a radical ideal

julia> Zred # now that we know the vanishing ideal, it is used for printing
Vanishing locus
  in Affine 2-space over Rational Field
  of ideal(y, x)

julia> Z = intersect(X,Y) # a non reduced scheme
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - y, y)

julia> Zred == Z
false

julia> Zred == reduced_scheme(Z)[1]
true

```
"""
function set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)
  Z = intersect(overlying_scheme(X), overlying_scheme(Y))
  # not sure how reduced vs geometrically reduced behaves hence check=true
  return AffineAlgebraicSet(Z, is_reduced=false, check=true)
end

########################################################
# (3) Closure of algebraic sets
########################################################

@doc raw"""
    closure(X::AbsAffineAlgebraicSet)

Return the closure of ``X`` in its ambient affine space.
"""
function closure(X::AbsAffineAlgebraicSet)
  Xcl = closure(overlying_scheme(X), ambient_space(X))
  return affine_algebraic_set(Xcl, check=false)
end

########################################################
# (4) Irreducible Components
########################################################

@doc raw"""
    irreducible_components(X::AbsAffineAlgebraicSet) -> Vector{AffineVariety}

Return the irreducible components of ``X`` defined over the base field of ``X``.

Note that they may be reducible over the algebraic closure.
See also [`geometric_irreducible_components`](@ref).
"""
function irreducible_components(X::AbsAffineAlgebraicSet)
  error("not implemented")
end

# special case for affine space
irreducible_components(X::AbsAffineAlgebraicSet{<:Field, MPolyRing}) = [X]

function irreducible_components(X::T) where {T <: AbsAffineAlgebraicSet{<:Field, <:MPolyQuoRing}}
  I = ideal(X)
  if is_one(I) # catch the empty set
    return T[]
  end
  P = minimal_primes(I)
  return T[vanishing_locus(p, is_reduced=true, check=false) for p in P]
end

@doc raw"""
    geometric_irreducible_components(X::AbsAffineAlgebraicSet)

Return the geometrically irreducible components of ``X``.

They are the irreducible components ``V_{ij}`` of ``X`` seen over an algebraically
closed field and given as a vector of tuples ``(A_i, V_{ij}, d_{ij})``, say,
where ``A_i`` is an algebraic set which is irreducible over the base field of ``X``
and ``V_{ij}`` represents a corresponding class of galois conjugated geometrically
irreducible components of ``A_i`` defined over a number field of degree
``d_{ij}`` whose generator prints as `_a`.

This is expensive and involves taking field extensions.
"""
function geometric_irreducible_components(X::AbsAffineAlgebraicSet)
  error("not implemented")
end

geometric_irreducible_components(X::AbsAffineAlgebraicSet{<:Field, <:MPolyRing}) = [X]

function geometric_irreducible_components(X::AbsAffineAlgebraicSet{QQField, <:MPolyQuoRing})
  I = ideal(X)
  if is_one(I) # catch the empty set
    return []
  end
  Pabs = absolute_primary_decomposition(I)
  return [(vanishing_locus(p[2], is_reduced=true, check=false), affine_variety(p[3], check=false),p[4]) for p in Pabs]
end
