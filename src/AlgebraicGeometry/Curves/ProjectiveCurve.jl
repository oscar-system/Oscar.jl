
export ProjectiveCurve

@doc raw"""
    ProjectiveCurve(I::MPolyIdeal)

Given a homogeneous ideal `I` of Krull dimension 2, return the projective curve defined by `I`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> M = matrix(R, 2, 3, [w x y; x y z])
[w   x   y]
[x   y   z]

julia> V = minors(M, 2)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 w*y - x^2
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> TC = ProjectiveCurve(I)
Projective curve defined by the ideal(w*y - x^2, w*z - x*y, x*z - y^2)
```
"""
@attributes mutable struct ProjectiveCurve{BaseRingType<:Field, RingType<:Ring} <: AbsProjectiveCurve{BaseRingType, RingType}
    X::ProjectiveAlgebraicSet{BaseRingType, RingType}

    function ProjectiveCurve(X::ProjectiveAlgebraicSet{S,T}; check::Bool=true) where {S,T}
        @check dim(X) == 1 || error("not of dimension one")
        new{S,T}(X)
    end
end

ProjectiveCurve(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true, is_radical::Bool=true) = ProjectiveCurve(algebraic_set(I;check=check,is_radical=is_radical);check=check)

projective_curve(I;kwargs...) = ProjectiveCurve(I;kwargs...)

underlying_scheme(X::ProjectiveCurve) = X.X
fat_scheme(X::ProjectiveCurve) = fat_scheme(underlying_scheme(X))

#=
function Base.show(io::IO, ::MIME"text/plain", C::ProjectivePlaneCurve)
  io = pretty(io)
  println(io, "Projective curve")
  print(io, Indent(), "defined by 0 = ", defining_equation(C), Dedent())
end
=#

################################################################################

@doc raw"""
    invert_birational_map(phi::Vector{T}, C::ProjectiveCurve) where {T <: MPolyRingElem}

Return a dictionary where `image` represents the image of the birational map
given by `phi`, and `inverse` represents its inverse, where `phi` is a
birational map of the projective curve `C` to its image in the projective
space of dimension `size(phi) - 1`.
Note that the entries of `inverse` should be considered as
representatives of elements in `R/image`, where `R` is the basering.
"""
function invert_birational_map(phi::Vector{T}, C::ProjectiveCurve) where {T <: MPolyRingElem}
    s = parent(phi[1])
    I = ideal(s, phi)
    singular_assure(I)
    IC = defining_ideal(C)
    singular_assure(IC)
    L = Singular.LibParaplanecurves.invertBirMap(I.gens.S, IC.gens.S)
    R = _fromsingular_ring(L[1])
    J = L[2][:J]
    psi = L[2][:psi]
    return Dict(:image => gens(ideal(R, J)), :inverse => gens(ideal(R, psi)))
end


################################################################################
# The union of two curves is a curve
Base.union(C::T, D::T) where T <: AbsProjectiveCurve = ProjectiveCurve(union(underlying_scheme(C),underlying_scheme(D)))
