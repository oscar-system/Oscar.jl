export ProjCurve, defining_ideal, curve_components, reduction, is_irreducible,
       jacobi_ideal

import Oscar.defining_ideal
################################################################################
abstract type ProjectiveCurve end
################################################################################
@doc Markdown.doc"""
    ProjCurve(I::MPolyIdeal)

Given a homogeneous ideal `I` of Krull dimension 2, return the projective curve defined by `I`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> M = matrix(R, 2, 3, [w x y; x y z])
[w   x   y]
[x   y   z]

julia> V = minors(M, 2)
3-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 w*y - x^2
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> TC = ProjCurve(I)
Projective curve defined by the ideal(w*y - x^2, w*z - x*y, x*z - y^2)
```
"""
mutable struct ProjCurve <: ProjectiveCurve
    I::MPolyIdeal
    ambiant_dim::Int
    components::Dict{ProjCurve, ProjCurve}
    function ProjCurve(I::MPolyIdeal)
        dim(I) == 2 || error("wrong dimension for a projective curve")
        n = nvars(base_ring(I)) - 1
        new(I, n, Dict{ProjCurve, ProjCurve}())
    end   
    function Base.show(io::IO, C::ProjCurve)
        if !get(io, :compact, false)
            print(io, "Projective curve defined by the ", C.I)
        else
            print(io, C.I)
        end
    end
end

################################################################################

@doc Markdown.doc"""
    defining_ideal(C::ProjCurve)

Return the defining ideal of the projective curve `C`.
"""
function defining_ideal(C::ProjCurve)
    return C.I
end

################################################################################

function Base.hash(C::ProjCurve, h::UInt)
  I = defining_ideal(C)
  return hash(I, h)
end

################################################################################
@doc Markdown.doc"""
    in(P::Oscar.Geometry.ProjSpcElem, C::ProjCurve)

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.

# Examples
```jldoctest
julia> S, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"])
(Multivariate Polynomial Ring in x, y, z, t over Rational Field, fmpq_mpoly[x, y, z, t])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z, t over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1]
  t -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z, t])

julia> I = ideal(T, [x^2, y^2*z, z^2])
ideal(x^2, y^2*z, z^2)

julia> C = Oscar.ProjCurve(I)
Projective curve defined by the ideal(x^2, y^2*z, z^2)


julia> PP = proj_space(QQ, 3)
(Projective space of dim 3 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2], x[3]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(2), QQ(0), QQ(5)])
(0 : 2 : 0 : 5)

julia> P in C
true
```
"""
function Base.in(P::Oscar.Geometry.ProjSpcElem, C::ProjCurve)
    I = defining_ideal(C)
    V = gens(I)
    return all(f -> iszero(evaluate(f, P.v)), V)
end

################################################################################
@doc Markdown.doc"""
    curve_components(C::ProjCurve)

Return a dictionary containing the irreducible components of `C` and the
corresponding reduced curve.
"""
function curve_components(C::ProjCurve)
    if isempty(C.components)
        I = defining_ideal(C)
        L = primary_decomposition(I)
        D = Dict(ProjCurve(q) => ProjCurve(p) for (q, p) in L if dim(q) == 2)
        C.components = D
    end
    return C.components
end

################################################################################

@doc Markdown.doc"""
    is_irreducible(C::ProjCurve)

Return `true` if `C` is irreducible, and `false` otherwise.

# Examples
```jldoctest
julia> S, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"])
(Multivariate Polynomial Ring in x, y, z, t over Rational Field, fmpq_mpoly[x, y, z, t])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z, t over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1]
  t -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z, t])

julia> I = ideal(T, [x^2, y^2*z, z^2])
ideal(x^2, y^2*z, z^2)

julia> C = Oscar.ProjCurve(I)
Projective curve defined by the ideal(x^2, y^2*z, z^2)

julia> Oscar.is_irreducible(C)
true
```
"""
function Oscar.is_irreducible(C::ProjCurve)
    L = curve_components(C)
    return length(L) == 1
end

################################################################################

@doc Markdown.doc"""
    reduction(C::ProjCurve)

Return the projective curve defined by the radical of the defining ideal of `C`.

# Examples
```jldoctest
julia> S, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"])
(Multivariate Polynomial Ring in x, y, z, t over Rational Field, fmpq_mpoly[x, y, z, t])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z, t over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1]
  t -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z, t])

julia> I = ideal(T, [x^2, y^2*z, z^2])
ideal(x^2, y^2*z, z^2)

julia> C = Oscar.ProjCurve(I)
Projective curve defined by the ideal(x^2, y^2*z, z^2)

julia> Oscar.reduction(C)
Projective curve defined by the ideal(z, x)
```
"""
function reduction(C::ProjCurve)
    J = radical(defining_ideal(C))
    return ProjCurve(J)
end

################################################################################

@doc Markdown.doc"""
    jacobi_ideal(C::ProjCurve)

Return the Jacobian ideal of the defining ideal of `C`.

# Examples
```jldoctest
julia> S, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"])
(Multivariate Polynomial Ring in x, y, z, t over Rational Field, fmpq_mpoly[x, y, z, t])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z, t over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1]
  t -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z, t])

julia> I = ideal(T, [x^2, y^2*z, z^2])
ideal(x^2, y^2*z, z^2)

julia> C = Oscar.ProjCurve(I)
Projective curve defined by the ideal(x^2, y^2*z, z^2)

julia> Oscar.jacobi_ideal(C)
ideal(4*x*y*z, 2*x*y^2, 4*x*z, 4*y*z^2)
```
"""
function Oscar.jacobi_ideal(C::ProjCurve)
    I = defining_ideal(C)
    R = base_ring(I)
    k = C.ambiant_dim - 1
    M = jacobi_matrix(gens(I))
    V = minors(M, k)
    filter!(!iszero, V)
    return ideal(R, V)
end

################################################################################

@doc Markdown.doc"""
    invert_birational_map(phi::Vector{T}, C::ProjCurve) where {T <: MPolyElem}

Return a dictionary where `image` represents the image of the birational map
given by `phi`, and `inverse` represents its inverse, where `phi` is a
birational map of the projective curve `C` to its image in the projective
space of dimension `size(phi) - 1`.
Note that the entries of `inverse` should be considered as
representatives of elements in `R/image`, where `R` is the basering.
"""
function invert_birational_map(phi::Vector{T}, C::ProjCurve) where {T <: MPolyElem}
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
