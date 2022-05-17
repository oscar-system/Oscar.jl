import Oscar: is_effective
export AffineCurveDivisor, ProjCurveDivisor, is_effective, multiplicity,
       divisor, curve, degree, divisor_ideals, global_sections,
       dimension_global_sections, is_linearly_equivalent, is_principal,
       curve_zero_divisor, principal_divisor

################################################################################

abstract type CurveDivisor end

################################################################################
# The curve C is assume to be smooth and irreducible, but it is not checked.
@doc Markdown.doc"""
    AffineCurveDivisor(C::AffinePlaneCurve{S}, D::Dict{Point{S}, Int}) where S <: FieldElem

Given a curve `C` which is assumed to be smooth and irreducible, return the divisor on the curve `C` defined by `D`.

# Examples
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^2 + y + x^2)
Affine plane curve defined by x^2 + y^2 + y

julia> P = Oscar.Point([QQ(0), QQ(0)])
Point with coordinates fmpq[0, 0]

julia> Q = Oscar.Point([QQ(0), QQ(-1)])
Point with coordinates fmpq[0, -1]

julia> Oscar.AffineCurveDivisor(C, Dict(P => 3, Q => -2))
3*fmpq[0, 0] - 2*fmpq[0, -1]
```
"""
struct AffineCurveDivisor{S <: FieldElem} <: CurveDivisor
    C::AffinePlaneCurve{S}
    divisor::Dict{Point{S}, Int}
    degree::Int
    function AffineCurveDivisor{S}(C::AffinePlaneCurve{S}, D::Dict{Point{S}, Int}) where S <: FieldElem
        for (P, m) in D
            P in C || error("The point ", P.coord, " is not on the curve")
            if m == 0
                delete!(D, P)
            end
        end
        new{S}(C, D, sum(values(D)))
    end
end

function AffineCurveDivisor(C::AffinePlaneCurve{S}) where S <: FieldElem
    return AffineCurveDivisor{S}(C, Dict{Point{S}, Int}())
end

function AffineCurveDivisor(C::AffinePlaneCurve{S}, D::Dict{Point{S}, Int}) where S <: FieldElem
    return AffineCurveDivisor{S}(C, D)
end

function AffineCurveDivisor(C::AffinePlaneCurve{S}, P::Point{S}, m::Int=1) where S <: FieldElem
    return AffineCurveDivisor{S}(C, Dict(P => m))
end

################################################################################

@doc Markdown.doc"""
    ProjCurveDivisor(C::ProjPlaneCurve{S}, D::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}) where S <: FieldElem

Given a curve `C` which is assumed to be smooth and irreducible, return the divisor on the curve `C` defined by `D`.

# Examples
```jldoctest
julia> S, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2 + y*z + x^2))
Projective plane curve defined by x^2 + y^2 + y*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
(0 : 0 : 1)

julia> Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-1), QQ(1)])
(0 : -1 : 1)

julia> D = Oscar.ProjCurveDivisor(C, Dict(P => 3, Q => -2))
3*(0 : 0 : 1) - 2*(0 : 1 : -1)
```
"""
mutable struct ProjCurveDivisor{S <: FieldElem} <: CurveDivisor
    C::ProjPlaneCurve{S}
    divisor::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}
    degree::Int
    prin_div::AbstractAlgebra.Generic.Frac{T} where T <: MPolyElem{S}
    function ProjCurveDivisor{S}(C::ProjPlaneCurve{S}, D::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}) where S <: FieldElem
        for (P, m) in D
            P in C || error("The point ", P, " is not on the curve")
            if m == 0
                delete!(D, P)
            end
        end
        new{S}(C, D, sum(values(D)))
    end
end

function ProjCurveDivisor(C::ProjPlaneCurve{S}) where S <: FieldElem
    return ProjCurveDivisor{S}(C, Dict{Oscar.Geometry.ProjSpcElem{S}, Int}())
end

function ProjCurveDivisor(C::ProjPlaneCurve{S}, D::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}) where S <: FieldElem
    return ProjCurveDivisor{S}(C, D)
end

function ProjCurveDivisor(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}, m::Int=1) where S <: FieldElem
    return ProjCurveDivisor{S}(C, Dict(P => m))
end

################################################################################
# To get a nice output for the divisors.

function AbstractAlgebra.expressify(@nospecialize(a::AffineCurveDivisor); context = nothing)
   prod = Expr(:call, :+)
   for (P, m) in a.divisor
      ep = AbstractAlgebra.expressify(P.coord; context = context)::String
      push!(prod.args, Expr(:call, :*, m, ep))
   end
   return prod
end

function AbstractAlgebra.expressify(@nospecialize(a::ProjCurveDivisor); context = nothing)
   prod = Expr(:call, :+)
   for (P, m) in a.divisor
      ep = AbstractAlgebra.expressify(P; context = context)
      push!(prod.args, Expr(:call, :*, m, ep))
   end
   return prod
end

function Base.show(io::IO, ::MIME"text/plain", a::CurveDivisor)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function Base.show(io::IO, a::CurveDivisor)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

################################################################################
################################################################################
@doc Markdown.doc"""
    curve(D::CurveDivisor)

Return the curve on which the divisor is considered.
"""
function curve(D::CurveDivisor)
    return D.C
end

@doc Markdown.doc"""
    degree(D::CurveDivisor)

Return the degree of the divisor.
"""
function Oscar.degree(D::CurveDivisor)
    return D.degree
end

################################################################################
@doc Markdown.doc"""
    curve_zero_divisor(C::ProjPlaneCurve{S}) where S <: FieldElem

Return the divisor `0` on the curve `C`.
"""
function curve_zero_divisor(C::ProjPlaneCurve{S}) where S <: FieldElem
    return ProjCurveDivisor{S}(C, Dict{Oscar.Geometry.ProjSpcElem{S}, Int}())
end

@doc Markdown.doc"""
    curve_zero_divisor(C::AffinePlaneCurve{S}) where S <: FieldElem

Return the divisor `0` on the curve `C`.
"""
function curve_zero_divisor(C::AffinePlaneCurve{S}) where S <: FieldElem
    return AffineCurveDivisor{S}(C, Dict{Point{S}, Int}())
end

################################################################################

function _check_same_curve(D::CurveDivisor, E::CurveDivisor)
    D.C == E.C || error("The divisors are not on the same curve")
end

################################################################################
################################################################################
# Sum of two divisors

function Base.:+(D::T, E::T) where T <: CurveDivisor
    _check_same_curve(D, E)
    return T(D.C, merge(+, D.divisor, E.divisor))
end

################################################################################

function _help_minus(D::CurveDivisor, E::CurveDivisor)
    _check_same_curve(D, E)
    F = copy(D.divisor)
    for (P, m) in E.divisor
        if haskey(F, P)
            F[P] -= m
        else
            F[P] = -m
        end
    end
    return F
end

function Base.:-(D::T, E::T) where T <: CurveDivisor
    return T(D.C, _help_minus(D, E))
end

################################################################################
################################################################################

function Base.:*(k::Int, D::T) where T <: CurveDivisor
    k == 0 && return curve_zero_divisor(D.C)
    k == 1 && return D
    return T(D.C, Dict((P, k*m) for (P, m) in D.divisor))
end

################################################################################

function ==(D::CurveDivisor, E::CurveDivisor)
    return D.C == E.C && D.divisor == E.divisor
end

################################################################################
################################################################################

@doc Markdown.doc"""
    is_effective(D::CurveDivisor)

Return `true` if `D` is an effective divisor, `false` otherwise.

# Examples
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^2 + y + x^2)
Affine plane curve defined by x^2 + y^2 + y

julia> P = Oscar.Point([QQ(0), QQ(0)])
Point with coordinates fmpq[0, 0]

julia> Q = Oscar.Point([QQ(0), QQ(-1)])
Point with coordinates fmpq[0, -1]

julia> D = Oscar.AffineCurveDivisor(C, Dict(P => 3, Q => -2))
3*fmpq[0, 0] - 2*fmpq[0, -1]

julia> Oscar.is_effective(D)
false
```
"""
function is_effective(D::CurveDivisor)
    return all(v >= 0 for v in values(D.divisor))
end

################################################################################
################################################################################
# TODO: For multiplicity and divisor, when it will be possible, take as an input
# an element of the (fraction field of the ) coordinate ring of the curve, and
# remove the curve from the input. It would be possible for the moment for
# polynomials but not for fractions. For homogeneity, we keep for the moment
# the following presentation.
###############################################################################

@doc Markdown.doc"""
    multiplicity(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Point{S}) where {S <: FieldElem, T <: MPolyElem{S}}
    multiplicity(C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Oscar.Geometry.ProjSpcElem{S})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}

Return the multiplicity of the rational function `phi` on the curve `C` at the point `P`.

# Examples
```jldoctest
julia> S, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2 + y*z + x^2))
Projective plane curve defined by x^2 + y^2 + y*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
(0 : 0 : 1)

julia> phi = T(x)//T(y)
x//y

julia> Oscar.multiplicity(C, phi, P)
-1
```
"""
function multiplicity(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Point{S}) where {S <: FieldElem, T <: MPolyElem{S}}
    f = divrem(phi.num, C.eq)
    g = divrem(phi.den, C.eq)
    !iszero(g[2]) || error("This is not a rational function on `C`")
    return multiplicity(C, f[2], P) - multiplicity(C, g[2], P)
end

function multiplicity(C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Oscar.Geometry.ProjSpcElem{S})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
    g = divrem(phi.den.f, C.eq.f)
    !iszero(g[2]) || error("This is not a rational function on the curve")
    f = divrem(phi.num.f, C.eq.f)
    !iszero(f[2]) || error("The numerator is zero on the curve")
    R = parent(C.eq)
    return multiplicity(C, R(f[2]), P) - multiplicity(C, R(g[2]), P)
end

###############################################################################

@doc Markdown.doc"""
    multiplicity(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}, P::Point{S}) where S <: FieldElem
    multiplicity(C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the multiplicity of the polynomial `F` on the curve `C` at the point `P`.
"""
function multiplicity(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}, P::Point{S}) where S <: FieldElem
    f = divrem(F, C.eq)
    if iszero(f[2])
        return Inf
    elseif is_constant(f[2])
        return 0
    else
        return intersection_multiplicity(C, AffinePlaneCurve(f[2]), P)
    end
end

function multiplicity(C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
    f = divrem(F.f, defining_equation(C))
    if iszero(f[2])
        return Inf
    elseif is_constant(f[2])
        return 0
    else
        R = parent(C.eq)
        return intersection_multiplicity(C, ProjPlaneCurve(R(f[2])), P)
    end
end

################################################################################
################################################################################

@doc Markdown.doc"""
    divisor(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}) where S <: FieldElem

Return the divisor defined by the polynomial `F` on the curve `C`.
"""
function divisor(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}) where S <: FieldElem
    f = divrem(F, C.eq)
    !iszero(f[2]) || error("The polynomial is zero on the curve")
    E = Dict{Point{S}, Int}()
    if !is_constant(f[2])
        CC = AffinePlaneCurve(f[2])
        L = curve_intersect(C, CC)
        for P in L[2]
            E[P] = intersection_multiplicity(C, CC, P)
        end
    end
    return AffineCurveDivisor(C, E)
end

################################################################################

@doc Markdown.doc"""
    divisor(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}) where {S <: FieldElem, T <: MPolyElem{S}}

Return the divisor defined by the rational function `phi` on the curve `C`.
"""
function divisor(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}) where {S <: FieldElem, T <: MPolyElem{S}}
    g = divrem(phi.den, C.eq)
    !iszero(g[2]) || error("This is not a rational function on the curve")
    f = divrem(phi.num, C.eq)
    !iszero(f[2]) || error("The numerator is zero on the curve")
    return divisor(C, f[2]) - divisor(C, g[2])
end

################################################################################

@doc Markdown.doc"""
    divisor([PP::Oscar.Geometry.ProjSpc{S}], C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}) where S <: FieldElem

Return the divisor defined by the polynomial `F` on the curve `C`. The points of the divisor are in the projective space `PP` if specified, or in a new projective space otherwise.
"""
function divisor(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}) where S <: FieldElem
    f = divrem(F.f, defining_equation(C))
    !iszero(f[2]) || error("The polynomial is zero on the curve")
    E = Dict{Oscar.Geometry.ProjSpcElem{S}, Int}()
    if !is_constant(f[2])
        R = parent(C.eq)
        CC = ProjPlaneCurve(R(f[2]))
        L = curve_intersect(PP, C, CC)
        for P in L[2]
            E[P] = intersection_multiplicity(C, CC, P)
        end
    end
    return ProjCurveDivisor(C, E)
end

function divisor(C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}) where S <: FieldElem
    R = parent(C.eq)
    PP = proj_space(R.R.base_ring, 2)
    return divisor(PP[1], C, F)
end

################################################################################

@doc Markdown.doc"""
    divisor(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}

Return the divisor defined by the rational function `phi` on the curve `C`.

# Examples
```jldoctest
julia> S, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2 + y*z + x^2))
Projective plane curve defined by x^2 + y^2 + y*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> phi = T(x)//T(y)
x//y

julia> Oscar.divisor(PP[1], C, phi)
-(0 : 0 : 1) + (0 : 1 : -1)
```
"""
function divisor(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
    g = divrem(phi.den.f, C.eq.f)
    !iszero(g[2]) || error("This is not a rational function on the curve")
    f = divrem(phi.num.f, C.eq.f)
    !iszero(f[2]) || error("The numerator is zero on the curve")
    R = parent(C.eq)
    return divisor(PP, C, R(f[2])) - divisor(PP, C, R(g[2]))
end

function divisor(C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
    R = parent(C.eq)
    PP = proj_space(R.R.base_ring, 2)
    return divisor(PP[1], C, phi)
end

################################################################################

################################################################################
# Ideal associated to a divisor

function divisor_ideals(D::ProjCurveDivisor)
    R = parent(defining_equation(D.C))
    V = gens(R)
    I = ideal(R, [R(1)])
    J = ideal(R, [R(1)])
    for (P, m) in D.divisor
        IP = ideal(R, [P.v[1]*V[2]-P.v[2]*V[1], P.v[1]*V[3]-P.v[3]*V[1], P.v[2]*V[3]-P.v[3]*V[2]])
        if m > 0
            I = intersect(I, IP^m)
        else
            J = intersect(J, IP^(-m))
        end
    end
    return [I, J]
end

################################################################################

################################################################################
################################################################################
# This code is copied from Singulars "divisors.lib", based on the Macaulay 2
# tutorial.

function jet(I::Oscar.MPolyIdeal, d::Int)
   singular_assure(I)
   R = base_ring(I)
   Is = Singular.jet(I.gens.S, d)
   return ideal(R, Is)
end

function _remove_zeros(I::T) where T <: Union{MPolyIdeal, MPolyQuoIdeal}
   R = base_ring(I)
   g = gens(I)
   filter!(x->x!= R(0), g)
   length(g) == 0 && return ideal(R, [R(0)])
   return ideal(R, g)
end

################################################################################
function _minimal_generating_set(I::Oscar.MPolyIdeal)
    R = base_ring(I)
    singular_assure(I)
    Ic = deepcopy(I.gens.S)
    Sx = I.gens.Sx
    Imin = Singular.Ideal(Sx, Singular.libSingular.idMinBase(Ic.ptr, Sx.ptr))
    return R.(gens(Imin))
end

# Remove components which are not codim 1
function _purify1(I::T, Q) where T <: Union{MPolyIdeal, MPolyQuoIdeal}
   R = base_ring(I)
   Id = _remove_zeros(I)
   gens(Id)[1] != R(0) || error("ideal assumed to be non-zero")
   IdQ = ideal(Q, gens(Id))
   f = ideal(Q, [gens(Id)[1]])
   res = f:(f:IdQ)
   Oscar.oscar_assure(res)
   return ideal(Q, _minimal_generating_set(res.I))
end

function _basis(I::Oscar.MPolyIdeal, d::Int)
   R = base_ring(I)
   m = ideal(R, gens(R))
   return _minimal_generating_set(_remove_zeros(jet(intersect(I, m^d),d)))
end

function _global_sections_helper(I::Oscar.MPolyIdeal, J::Oscar.MPolyIdeal, q::Oscar.MPolyIdeal)
   R = base_ring(I)
   R == base_ring(J) || error("base rings do not match")
   R == base_ring(q) || error("base rings do not match")
   Q, _ = quo(R, q)
   Ip = _purify1(I, Q)
   Jp = _purify1(J, Q)
   Is = _remove_zeros(Ip)
   f = gens(Is)[1]
   fJ = ideal(Q, [f*g for g in gens(Jp)])
   q1 = fJ:Ip
   P = _purify1(q1, Q)
   Oscar.oscar_assure(P)
   B = _basis(P.I, total_degree(f.f))
   arr = normal_form(B, q)
   LD = ideal(R, arr)
   return [_remove_zeros(LD), f.f]
end

################################################################################

function _global_sections_ideals(D::ProjCurveDivisor)
    F = defining_equation(D.C)
    R = parent(F)
    q = ideal(R, [F])
    E = divisor_ideals(D)
    _global_sections_helper(E[1], E[2], q)
end

################################################################################

@doc Markdown.doc"""
    global_sections(D::ProjCurveDivisor)

Return a set of generators of the global sections of the sheaf associated to the divisor `D` of a smooth and irreducible projective curve.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2*z - x*(x-z)*(x+3*z)))
Projective plane curve defined by -x^3 - 2*x^2*z + 3*x*z^2 + y^2*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
(0 : 1 : 0)

julia> D = Oscar.ProjCurveDivisor(C, P, 4)
4*(0 : 1 : 0)

julia> Oscar.global_sections(D)
4-element Vector{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}:
 1
 y//z
 x//z
 x^2//z^2
```
"""
function global_sections(D::ProjCurveDivisor)
    L = _global_sections_ideals(D)
    N = gens(L[1])
    return [g//L[2] for g in N]
end

################################################################################
@doc Markdown.doc"""
    dimension_global_sections(D::ProjCurveDivisor)

Return the dimension of the global sections of the sheaf associated to the divisor `D` of a smooth and irreducible projective curve.
"""
function dimension_global_sections(D::ProjCurveDivisor)
    L = global_sections(D)
    if length(L) != 1
        return length(L)
    elseif iszero(L[1])
        return 0
    else
        return 1
    end
end

################################################################################

function _section_ideal(f::Oscar.MPolyElem, g::Oscar.MPolyElem, D::ProjCurveDivisor)
    E = divisor_ideals(D)
    R = parent(f)
    F = defining_equation(D.C)
    R == parent(g) && R == parent(F) || error("base rings do not match")

    q = ideal(R, [F])
    Q, _ = quo(R, q)
    J = ideal(Q, [g])
    fE1 = ideal(Q, gens(f*E[1]))
    E2 = ideal(Q, gens(E[2]))
    return _purify1((fE1:J):E2, Q)
end

################################################################################

function _linearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
    F = D - E
    T = parent(D.C.eq)
    L = _global_sections_ideals(F)
    if length(gens(L[1])) !=1
        return T(0)//T(1)
    elseif iszero(gens(L[1])[1])
        return T(0)//T(1)
    else
        V = _section_ideal(gens(L[1])[1], L[2], F)
        Q = base_ring(V)
        if V == ideal(Q, [Q(1)])
            return T(gens(L[1])[1])//T(L[2])
        else
            return T(0)//T(1)
        end
    end
end

################################################################################
@doc Markdown.doc"""
    is_linearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)

Return `true` if the divisors `D` and `E` are linearly equivalent, and `false`
otherwise

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2*z - x*(x-z)*(x+3*z)))
Projective plane curve defined by -x^3 - 2*x^2*z + 3*x*z^2 + y^2*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
(0 : 1 : 0)

julia> R = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
(0 : 0 : 1)

julia> E = Oscar.ProjCurveDivisor(C, P)
(0 : 1 : 0)

julia> F = Oscar.ProjCurveDivisor(C, R)
(0 : 0 : 1)

julia> Oscar.is_linearly_equivalent(E, F)
false
```
"""
function is_linearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
    return !iszero(_linearly_equivalent(D, E))
end

################################################################################
@doc Markdown.doc"""
    is_principal(D::ProjCurveDivisor{S}) where S <: FieldElem

Return `true` if the divisor `D` is principal, and `false` otherwise

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2*z - x*(x-z)*(x+3*z)))
Projective plane curve defined by -x^3 - 2*x^2*z + 3*x*z^2 + y^2*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
(0 : 1 : 0)

julia> E = Oscar.ProjCurveDivisor(C, P)
(0 : 1 : 0)

julia> Oscar.is_principal(E)
false
```
"""
function Oscar.is_principal(D::ProjCurveDivisor{S}) where S <: FieldElem
    if !isdefined(D, :prin_div)
        D.prin_div = _linearly_equivalent(D, curve_zero_divisor(D.C))
    end
    return !iszero(D.prin_div)
end

################################################################################
@doc Markdown.doc"""
    principal_divisor(D::ProjCurveDivisor{S}) where S <: FieldElem

If the divisor `D` is principal, return a rational function `phi` such that `D`
is linearly equivalent to the divisor defined by `phi`.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2*z - x*(x-z)*(x+3*z)))
Projective plane curve defined by -x^3 - 2*x^2*z + 3*x*z^2 + y^2*z

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
(0 : 1 : 0)

julia> R = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
(0 : 0 : 1)

julia> E = Oscar.ProjCurveDivisor(C, P, 2)
2*(0 : 1 : 0)

julia> F = Oscar.ProjCurveDivisor(C, R, 2)
2*(0 : 0 : 1)

julia> G = 2*E - 2*F
-4*(0 : 0 : 1) + 4*(0 : 1 : 0)

julia> Oscar.principal_divisor(G)
x^2//z^2
```
"""
function principal_divisor(D::ProjCurveDivisor{S}) where S <: FieldElem
    is_principal(D) || error("the divisor is not principal")
    return D.prin_div
end
################################################################################
