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

This computes the radical of ``I`` if `check=true`.
Otherwise, Oscar takes on faith that ``I`` is radical.

```
julia> R, (x,y) = GF(2)[:x,:y];

julia> X = Oscar.vanishing_locus(ideal([y^2+y+x^3+1,x]))
Vanishing locus
  in Affine 2-space over GF(2)
  of ideal(x, y^2 + y + 1)

```
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

affine_algebraic_set(I::MPolyIdeal{<:MPolyElem}; check::Bool=true) = vanishing_locus(I, check=check)

@doc raw"""
    vanishing_locus(p::MPolyRingElem, check::Bool=true)

Return the vanishing locus of the multivariate polynomial `p`.

This computes the radical of ``I`` if `check=true`.
Otherwise Oscar takes on faith that ``I`` is radical.

```
julia> R, (x,y) = QQ[:x,:y];

julia> X = Oscar.vanishing_locus((y^2+y+x^3+1)*x^2)
Vanishing locus
  in Affine 2-space over QQ
  of ideal(x^4 + x*y^2 + x*y + x)

julia> R, (x,y) = GF(2)[:x,:y];

julia> X = Oscar.vanishing_locus((y^2+y+x^3+1)*x^2)
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
irreducible_components(X::AbsAffineAlgebraicSet{<:Field,MPolyRing}) = [X]

function irreducible_components(X::T) where {T <: AbsAffineAlgebraicSet{<:Field, <:MPolyQuoRing}}
  I = vanishing_ideal(X)
  if is_one(I) # catch the empty set
    return T[]
  end
  P = minimal_primes(I)
  return T[vanishing_locus(p, check=false) for p in P]
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
  I = vanishing_ideal(X)
  if is_one(I) # catch the empty set
    return []
  end
  Pabs = absolute_primary_decomposition(I)
  return [(vanishing_locus(p[2], check=false), affine_variety(p[3], check=false),p[4]) for p in Pabs]
end
