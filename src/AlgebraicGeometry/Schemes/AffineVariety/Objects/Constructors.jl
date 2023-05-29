########################################################
# (1) Generic constructors
########################################################
@doc raw"""
    affine_variety(X::Spec; check::Bool=true) -> AffineVariety

Convert ``X`` to an affine variety.

If check is set, then compute the reduced scheme of `X` first.
"""
function affine_variety(X::Spec{<:Field}; check::Bool=true)
  check  ||  AffineVariety(X, check=check)
  Xred,_ = reduced_scheme(X)
  return AffineVariety(Xred, check=check)
end


@doc raw"""
    affine_variety(I::MPolyIdeal; check=true) -> AffineVariety

Return the affine variety defined by the prime ideal ``I``.

Since our varieties are irreducible, we check that ``I`` stays prime when
viewed over the algebraic closure. This is an expensive check that can be disabled.

```jldoctest
julia> R, (x,y) = QQ[:x,:y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> affine_variety(ideal([x,y]))
Affine variety
 in Affine 2-space over QQ
defined by ideal(x, y)

```
Over fields different from `QQ`, currently, we cannot check for irreducibility
over the algebraic closure. But if you know that the ideal in question defines
a variety, you can construct it by disabling the check.
```jldoctest
julia> R, (x,y) = GF(2)[:x,:y];

julia> affine_variety(x^3+y+1,check=false)
Affine variety
 in Affine 2-space over GF(2)
defined by ideal(x^3 + y + 1)

```
"""
affine_variety(I::MPolyIdeal; check=true) = AffineVariety(Spec(quo(base_ring(I),I)[1]), check=check)

@doc raw"""
    affine_variety(R::Ring; check=true)

Return the affine variety with coordinate ring `R`.

We require that ``R`` is a finitely generated algebra over a field ``k`` and
moreover that the base change of ``R`` to the algebraic closure ``\bar k``
is an integral domain.

```jldoctest
julia> R, (x,y) = QQ[:x,:y];

julia> Q,_ = quo(R,ideal([x,y]));

julia> affine_variety(Q)
Affine variety
 in Affine 2-space over QQ
defined by ideal(x, y)

```
"""
affine_variety(R::MPolyAnyRing; check=true) = AffineVariety(Spec(R), check=check)

@doc raw"""
    function affine_variety(f::MPolyRingElem{<:Field}; check::Bool=true)

Return the affine variety defined as the vanishing locus of the multivariate polynomial `f`.

This checks that `f` is irreducible over the algebraic closure.

```jldoctest
julia> A2 = affine_space(QQ,[:x,:y]);

julia> (x,y) = coordinates(A2);

julia> affine_variety(y^2-x^3-1)
Affine variety
 in Affine 2-space over QQ
defined by ideal(-x^3 + y^2 - 1)

```
"""
function affine_variety(f::MPolyRingElem{<:FieldElem}; check::Bool=true)
  if check
    is_irreducible(f) || error("polynomial is reducible")
    ff = factor_absolute(f)[2]
    # deal with weird type instability in factor_absolute
    @assert ff[2] == 1
    g = ff[1]
    (length(g) == 1 || (length(g)==2 && isone(g[2]))) || error("polynomial is not absolutely irreducible")
  end
  return affine_variety(ideal([f]), check=false)
end



########################################################
# (3) Closure of varieties in their ambient space are varieties
########################################################

function closure(X::AffineVariety)
  Xcl = closure(X, ambient_space(X))
  return affine_algebraic_set(Xcl, check=false)
end


