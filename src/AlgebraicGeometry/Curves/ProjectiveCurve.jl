

@doc raw"""
    ProjectiveCurve

A reduced projective curve, defined as the vanishing locus of
a homogeneous (but not necessarily radical) ideal.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> M = matrix(R, 2, 3, [w x y; x y z])
[w   x   y]
[x   y   z]

julia> V = minors(M, 2)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 w*y - x^2
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> TC = projective_curve(I)
Projective curve
  in projective 3-space over QQ with coordinates [w, x, y, z]
defined by ideal (w*y - x^2, w*z - x*y, x*z - y^2)

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


function Base.show(io::IO, ::MIME"text/plain", X::ProjectiveCurve)
  io = pretty(io)
  println(io, "Projective curve")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  if isdefined(X, :Xred)
    I = vanishing_ideal(X)
  else
    I = fat_ideal(X)
  end
  print(io, Dedent(), "defined by ", Lowercase(), I)
end

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
    IC = defining_ideal(C)
    L = Singular.LibParaplanecurves.invertBirMap(singular_generators(I), singular_generators(IC))
    R = _fromsingular_ring(L[1])
    J = L[2][:J]
    psi = L[2][:psi]
    return Dict(:image => gens(ideal(R, J)), :inverse => gens(ideal(R, psi)))
end


################################################################################
# The union of two curves is a curve
Base.union(C::T, D::T) where T <: AbsProjectiveCurve = ProjectiveCurve(union(underlying_scheme(C),underlying_scheme(D)))
