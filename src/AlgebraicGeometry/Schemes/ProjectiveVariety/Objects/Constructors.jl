########################################################
# (1) Generic constructors
########################################################
@doc raw"""
    projective_variety(X::AbsProjectiveScheme; check::Bool=true) -> ProjectiveVariety

Convert ``X`` to a projective variety.

If `check` is set, compute the reduced scheme of `X` first.
"""
function projective_variety(X::AbsProjectiveScheme{<:Field}; check::Bool=true)
  check  ||  return ProjectiveVariety(X, check=check)
  Xred = reduced_scheme(X)
  return ProjectiveVariety(Xred, check=check)
end

# documented above
function projective_variety(X::AbsProjectiveAlgebraicSet{<:Field}; check::Bool=true)
  check  ||  return ProjectiveVariety(X, check=check)
  Xred = X # already geometrically reduced
  return ProjectiveVariety(Xred, check=check)
end

@doc raw"""
    projective_variety([R::MPolyDecRing,] I::MPolyIdeal; check::Bool=true) -> ProjectiveVariety

Return the projective variety defined by the homogeneous prime ideal ``I``.

Since in our terminology varieties are irreducible over the algebraic closure,
we check that ``I`` stays prime when viewed over the algebraic closure.
This is an expensive check that can be disabled.
Note that the ideal ``I`` must live in a standard graded ring.

```jldoctest
julia> P3 = projective_space(QQ,3)
Projective space of dimension 3
  with homogeneous coordinates s0 s1 s2 s3
  over Rational Field

julia> (s0,s1,s2,s3) = homogeneous_coordinates(P3);

julia> X = projective_variety(s0^3 + s1^3 + s2^3 + s3^3)
Projective variety
  in Projective 3-space over Rational Field
  defined by ideal(s0^3 + s1^3 + s2^3 + s3^3)

julia> dim(X)
2

julia> Y = projective_variety(ideal([s0^3 + s1^3 + s2^3 + s3^3, s0]))
Projective variety
  in Projective 3-space over Rational Field
  defined by ideal(s0^3 + s1^3 + s2^3 + s3^3, s0)

julia> dim(Y)
1
```
"""
projective_variety(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true) = ProjectiveVariety(ProjectiveScheme(base_ring(I),I), check=check)

function projective_variety(R::MPolyDecRing, I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
  @req base_ring(I) === R "ideal must be defined over R"
  ProjectiveVariety(ProjectiveScheme(R,I), check=check)
end


@doc raw"""
    projective_variety(R::Ring; check::Bool=true)

Return the projective variety defined by the ``\mathbb{Z}`` standard graded ring ``R``.

We require that ``R`` is a finitely generated algebra over a field ``k`` and
moreover that the base change of ``R`` to the algebraic closure ``\bar k``
is an integral domain.
"""
projective_variety(R::Ring; check::Bool=true) = ProjectiveVariety(ProjectiveScheme(R), check=check)

projective_variety(R::MPolyDecRing; check::Bool=true) = ProjectiveVariety(ProjectiveScheme(R), check=check)


@doc raw"""
    projective_variety(f::MPolyDecRingElem; check=true)

Return the projective variety defined by the homogeneous polynomial `f`.

This checks that `f` is absolutely irreducible.
"""
function projective_variety(f::MPolyDecRingElem; check=true)
  if check
    is_irreducible(f) || error("polynomial is reducible")
    ff = factor_absolute(forget_decoration(f))[2]
    # deal with weird type instability in factor_absolute
    # https://github.com/oscar-system/Oscar.jl/issues/2211
    @assert ff[2] == 1
    g = ff[1]
    (length(g) == 1 || (length(g)==2 && isone(g[2]))) || error("polynomial is not absolutely irreducible")
  end
  return projective_variety(ideal([f]), check=false)
end


########################################################
# (1) projective space
########################################################

function projective_space(A::Field, var_symb::Vector{VarName})
  n = length(var_symb)
  R, _ = polynomial_ring(A, var_symb)
  S, _ = grade(R, [1 for i in 1:n ])
  return projective_variety(projective_scheme(S), check=false)
end


function projective_space(
    A::CoeffRingType,
    r::Int;
    var_name::VarName=:s
  ) where {CoeffRingType<:Field}
  R, _ = polynomial_ring(A, [Symbol(var_name,i) for i in 0:r])
  S, _ = grade(R, [1 for i in 0:r ])
  return projective_variety(projective_scheme(S), check=false)
end
