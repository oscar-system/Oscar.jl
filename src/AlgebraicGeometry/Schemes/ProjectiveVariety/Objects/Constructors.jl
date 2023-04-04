########################################################
# (1) Generic constructors
########################################################
@doc Markdown.doc"""
    projective_variety(X::AbsProjectiveScheme; check::Bool=true) -> ProjectiveVariety

Convert ``X`` to an projective variety.

If check is set, then compute the reduced scheme of `X` first.
"""
function projective_variety(X::AbsProjectiveScheme{<:Field}; check::Bool=true)
  check  ||  ProjectiveVariety(X, check=check)
  Xred,_ = reduced_scheme(X)
  return ProjectiveVariety(Xred, check=check)
end

@doc Markdown.doc"""
    projective_variety(I::MPolyIdeal; check::Bool=true) -> ProjectiveVariety

Return the projective variety defined by the homogeneous prime ideal ``I``.

Since our varieties are irreducible, we check that ``I`` stays prime when
viewed over the algebraic closure. This is an expensive check that can be disabled.
"""
projective_variety(I::MPolyIdeal; check::Bool=true) = ProjectiveVariety(ProjectiveScheme(base_ring(I),I), check=check)

@doc Markdown.doc"""
    projective_variety(R::Ring; check::Bool=true)

Return the projective variety defined by ``R``.

We require that ``R`` is a finitely generated algebra over a field ``k`` and
moreover that the base change of ``R`` to the algebraic closure ``\bar k``
is an integral domain.
"""
projective_variety(R::Ring; check::Bool=true) = ProjectiveVariety(ProjectiveScheme(R), check=check)


function projective_variety(f::MPolyDecRingElem; check=true)
  if check
    is_irreducible(f) || error("polynomial is reducible")
    ff = factor_absolute(f)[2]
    # deal with weird type instability in factor_absolute
    @assert ff[2] == 1
    g = ff[1]
    (length(g) == 1 || (length(g)==2 && isone(g[2]))) || error("polynomial is not absolutely irreducible")
  end
  return projective_variety(ideal([f]), check=false)
end


