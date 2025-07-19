########################################################
# (1) Generic constructors
########################################################
@doc raw"""
    variety(X::AbsAffineScheme; is_reduced::false, check::Bool=true) -> AffineVariety

Convert ``X`` to an affine variety.

If `is_reduced` is set, assume that `X` is already reduced.
"""
function variety(X::AbsAffineScheme{<:Field}; is_reduced=false, check::Bool=true)
  X = algebraic_set(X, is_reduced=is_reduced, check=check)
  return variety(X, check=check)
end

function variety(X::AbsAffineAlgebraicSet; check::Bool=true)
  return AffineVariety(X, check=check)
end

@doc raw"""
    variety(I::MPolyIdeal; check=true) -> AffineVariety

Return the affine variety defined by the ideal ``I``.

By our convention, varieties are absolutely irreducible.
Hence we check that the radical of ``I`` is prime and stays prime when
viewed over the algebraic closure. This is an expensive check that can be disabled.

```jldoctest
julia> R, (x,y) = QQ[:x,:y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> variety(ideal([x,y]))
Affine variety
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (x, y)

```
Over fields different from `QQ`, currently, we cannot check for irreducibility
over the algebraic closure. But if you know that the ideal in question defines
a variety, you can construct it by disabling the check.
```jldoctest
julia> R, (x,y) = GF(2)[:x,:y];

julia> variety(x^3+y+1, check=false)
Affine variety
  in affine 2-space over GF(2) with coordinates [x, y]
defined by ideal (x^3 + y + 1)

```
"""
function variety(I::MPolyIdeal; is_radical=false, check=true)
  X = algebraic_set(I, is_radical=is_radical, check=check)
  return AffineVariety(X ,check=check)
end


@doc raw"""
    variety(R::Ring; check=true)

Return the affine variety with coordinate ring `R`.

We require that ``R`` is a finitely generated algebra over a field ``k`` and
moreover that the base change of ``R`` to the algebraic closure
``\bar k`` is an integral domain.

```jldoctest
julia> R, (x,y) = QQ[:x,:y];

julia> Q,_ = quo(R,ideal([x,y]));

julia> variety(Q)
Affine variety
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (x, y)

```
"""
variety(R::MPolyAnyRing; check=true) = variety(spec(R), check=check)

@doc raw"""
    variety(f::MPolyRingElem{<:Field}; check::Bool=true)

Return the affine variety defined by the multivariate polynomial `f`.

This checks that `f` is irreducible over the algebraic closure.

```jldoctest
julia> A2 = affine_space(QQ,[:x,:y])
Affine space of dimension 2
  over rational field
with coordinates [x, y]

julia> (x,y) = coordinates(A2);

julia> variety(y^2-x^3-1)
Affine variety
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (x^3 - y^2 + 1)
```
"""
function variety(f::MPolyRingElem{<:FieldElem}; check::Bool=true)
  if check
    is_irreducible(f) || error("polynomial is reducible")
    ff = factor_absolute(f)[2]
    # deal with weird type instability in factor_absolute
    @assert ff[2] == 1
    g = ff[1]
    (length(g) == 1 || (length(g)==2 && isone(g[2]))) || error("polynomial is not absolutely irreducible")
  end
  return variety(ideal([f]), check=false)
end



########################################################
# (3) Closure of varieties in their ambient space are varieties
########################################################

function closure(X::AffineVariety)
  Xcl = closure(X, ambient_space(X))
  return algebraic_set(Xcl, check=false)
end


