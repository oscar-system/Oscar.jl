########################################################
# (1) Generic constructors
########################################################
@doc raw"""
    variety(X::AbsProjectiveScheme; is_reduced::Bool=false, check::Bool=true) -> ProjectiveVariety

Convert ``X`` to a projective variety by considering its reduced structure
"""
function variety(X::AbsProjectiveScheme{<:Field}; is_reduced::Bool=false, check::Bool=true)
  Xset = algebraic_set(X, is_reduced=is_reduced, check=check)
  return ProjectiveVariety(Xset, check=check)
end


@doc raw"""
    variety(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true, is_radical::Bool=false) -> ProjectiveVariety

Return the projective variety defined by the homogeneous prime ideal ``I``.

Since in our terminology varieties are irreducible over the algebraic closure,
we check that ``I`` stays prime when viewed over the algebraic closure.
This is an expensive check that can be disabled.
Note that the ideal ``I`` must live in a standard graded ring.

```jldoctest
julia> P3 = projective_space(QQ,3)
Projective space of dimension 3
  over rational field
with homogeneous coordinates [s0, s1, s2, s3]

julia> (s0,s1,s2,s3) = homogeneous_coordinates(P3);

julia> X = variety(s0^3 + s1^3 + s2^3 + s3^3)
Projective variety
  in projective 3-space over QQ with coordinates [s0, s1, s2, s3]
defined by ideal (s0^3 + s1^3 + s2^3 + s3^3)

julia> dim(X)
2

julia> Y = variety(ideal([s0^3 + s1^3 + s2^3 + s3^3, s0]))
Projective variety
  in projective 3-space over QQ with coordinates [s0, s1, s2, s3]
defined by ideal (s0, s1^3 + s2^3 + s3^3)

julia> dim(Y)
1
```
"""
function variety(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true,
                 is_radical::Bool=false)
    X = algebraic_set(I, check=check, is_radical=is_radical)
    return ProjectiveVariety(X, check=check)
end

@doc raw"""
    variety(R::GradedRing; check::Bool=true)

Return the projective variety defined by the ``\mathbb{Z}`` standard graded ring ``R``.

We require that ``R`` is a finitely generated algebra over a field ``k`` and
moreover that the base change of ``R`` to the algebraic closure ``\bar k``
is an integral domain.
"""
variety(R::Ring; check::Bool=true) = ProjectiveVariety(proj(R), check=check)

variety(R::MPolyDecRing; check::Bool=true) = ProjectiveVariety(proj(R), check=check)


@doc raw"""
    variety(f::MPolyDecRingElem; check=true)

Return the projective variety defined by the homogeneous polynomial `f`.

This checks that `f` is absolutely irreducible.
"""
function variety(f::MPolyDecRingElem; check::Bool=true, is_radical::Bool=false)
  @check begin
    is_irreducible(f) || error("polynomial is reducible")
    ff = factor_absolute(forget_decoration(f))[2]
    # deal with weird type instability in factor_absolute
    # https://github.com/oscar-system/Oscar.jl/issues/2211
    @assert ff[2] == 1
    g = ff[1]
    (length(g) == 1 || (length(g)==2 && isone(g[2]))) || error("polynomial is not absolutely irreducible")
    true
  end
  return variety(ideal([f]), check=false, is_radical=is_radical)
end


########################################################
# (1) projective space
########################################################

function projective_space(A::Field, var_symb::Vector{VarName})
  n = length(var_symb)
  S, _ = graded_polynomial_ring(A, var_symb; cached=false)
  return variety(proj(S), check=false)
end


function projective_space(
    A::CoeffRingType,
    r::Int;
    var_name::VarName=:s
  ) where {CoeffRingType<:Field}
  S, _ = graded_polynomial_ring(A, [Symbol(var_name,i) for i in 0:r]; cached=false)
  return variety(proj(S), check=false)
end
