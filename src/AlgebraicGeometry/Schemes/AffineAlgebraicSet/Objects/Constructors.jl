########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    algebraic_set(X::AffineScheme; is_reduced::Bool=false, check::Bool=true) -> AffineAlgebraicSet

Convert `X` to an `AffineAlgebraicSet` by considering its reduced structure.

If `is_reduced` is set, assume that `X` is already reduced.
If `is_reduced` and `check` are set,
check that `X` is actually geometrically reduced as claimed.
"""
function algebraic_set(X::AffineScheme; is_reduced::Bool=false, check::Bool=true)
  return AffineAlgebraicSet(X; is_reduced, check)
end

@doc raw"""
    algebraic_set(I::MPolyIdeal; is_radical::Bool=false, check::Bool=true)

Return the affine algebraic set defined by ``I``.

If `is_radical` is set, assume that ``I`` is a radical ideal.
```jldoctest
julia> R, (x,y) = GF(2)[:x,:y];

julia> X = algebraic_set(ideal([y^2+y+x^3+1,x]))
Affine algebraic set
  in affine 2-space over GF(2) with coordinates [x, y]
defined by ideal (x^3 + y^2 + y + 1, x)

```
"""
function algebraic_set(I::MPolyIdeal{<:MPolyRingElem}; is_radical::Bool=false, check::Bool=true)
  X = spec(base_ring(I), I)
  return algebraic_set(X; is_reduced=is_radical, check)
end

@doc raw"""
    algebraic_set(p::MPolyRingElem)

Return the affine algebraic set defined by the multivariate polynomial `p`.

```jldoctest
julia> R, (x,y) = QQ[:x,:y];

julia> X = algebraic_set((y^2+y+x^3+1)*x^2)
Affine algebraic set
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (x^5 + x^2*y^2 + x^2*y + x^2)

julia> R, (x,y) = GF(2)[:x,:y];

julia> X = algebraic_set((y^2+y+x^3+1)*x^2)
Affine algebraic set
  in affine 2-space over GF(2) with coordinates [x, y]
defined by ideal (x^5 + x^2*y^2 + x^2*y + x^2)

```
"""
function algebraic_set(p::MPolyRingElem; is_radical::Bool =false, check::Bool=true)
  I = ideal(parent(p), p)
  return algebraic_set(I; check, is_radical)
end

########################################################
# (2) Intersections of algebraic sets
########################################################

@doc raw"""
    set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)

Return the set theoretic intersection of `X` and `Y` as an algebraic set.

```jldoctest set_theoretic_intersection
julia> A = affine_space(QQ, [:x,:y])
Affine space of dimension 2
  over rational field
with coordinates [x, y]

julia> (x, y) = coordinates(A)
2-element Vector{QQMPolyRingElem}:
 x
 y

julia> X = algebraic_set(ideal([y - x^2]))
Affine algebraic set
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (-x^2 + y)

julia> Y = algebraic_set(ideal([y]))
Affine algebraic set
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (y)

julia> Zred = set_theoretic_intersection(X, Y)
Affine algebraic set
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (-x^2 + y, y)


```
Note that the set theoretic intersection forgets the intersection multiplicities
which the scheme theoretic intersection remembers. Therefore they are different.

```jldoctest set_theoretic_intersection
julia> Z = intersect(X, Y) # a non reduced scheme
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal (x^2 - y, y)

julia> Zred == Z
false

julia> Zred == reduced_scheme(Z)[1]
true

```
"""
function set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)
  Z = intersect(fat_scheme(X), fat_scheme(Y))
  return algebraic_set(Z)
end

function union(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)
  Z = union(fat_scheme(X), fat_scheme(Y))
  return algebraic_set(Z)
end

########################################################
# (3) Closure of algebraic sets
########################################################

@doc raw"""
    closure(X::AbsAffineAlgebraicSet)

Return the closure of ``X`` in its ambient affine space.
"""
function closure(X::AbsAffineAlgebraicSet)
  Xcl = closure(fat_scheme(X), ambient_space(X))
  return algebraic_set(Xcl, check=false)
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

function irreducible_components(X::AbsAffineAlgebraicSet{S,T}) where {S<:Field, T<:MPolyQuoRing}
  I = fat_ideal(X)
  P = minimal_primes(I)
  if length(P)==1 && is_one(P[1]) # catch the empty set for now :-(
    return AffineAlgebraicSet{S,T}[]
  end
  return AffineAlgebraicSet{S,T}[algebraic_set(p, is_radical=true, check=false) for p in P]
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

function geometric_irreducible_components(X::AbsAffineAlgebraicSet{QQField, T} where {T<:MPolyQuoRing})
  I = vanishing_ideal(X)
  if is_one(I) # catch the empty set ... not typestable anyways
    return AffineAlgebraicSet{QQField, T}[]
  end
  Pabs = absolute_primary_decomposition(I)
  return [(algebraic_set(p[2], is_radical=true, check=false), variety(p[3], check=false),p[4]) for p in Pabs]
end
