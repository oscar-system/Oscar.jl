########################################################
# (1) Generic constructors
########################################################
@doc Markdown.doc"""
    affine_variety(X::Spec; check::Bool=true) -> AffineVariety

Convert ``X`` to an affine variety.

If check is set, then compute the reduced scheme of `X` first.
"""
function affine_variety(X::Spec{<:Field}; check::Bool=true)
  check  ||  AffineVariety(X, check=check)
  Xred,_ = reduced_scheme(X)
  return AffineVariety(Xred, check=check)
end

@doc Markdown.doc"""
    affine_variety(I::MPolyIdeal; check=true) -> AffineVariety

Return the affine variety defined by the prime ideal ``I``.

Since our varieties are irreducible, we check that ``I`` stays prime when
viewed over the algebraic closure. This is an expensive check that can be disabled.
"""
affine_variety(I::MPolyIdeal; check=true) = AffineVariety(Spec(quo(base_ring(I),I)[1]), check=check)

@doc Markdown.doc"""
    affine_variety(R::Ring; check=true)

Return the affine variety defined by ``R``.

We require that ``R`` is a finitely generated algebra over a field ``k`` and
moreover that the base change of ``R`` to the algebraic closure ``\bar k``
is an integral domain.
"""
affine_variety(R::MPolyAnyRing; check=true) = AffineVariety(Spec(R), check=check)


function affine_variety(f::MPolyRingElem; check=true)
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
# (3) Closure of algebraic sets in their ambient space
########################################################

function closure(X::AffineVariety)
  Xcl = closure(X, ambient_space(X))
  return affine_algebraic_set(Xcl, check=false)
end


