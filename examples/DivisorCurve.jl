export AffineCurveDivisor, ProjCurveDivisor, iseffective, multiplicity,
       divisor, curve, degree

################################################################################

abstract type CurveDivisor end

################################################################################
# The curve C is assume to be smooth and irreducible, but it is not checked.
@doc Markdown.doc"""
    AffineCurveDivisor(C::AffinePlaneCurve{S}, D::Dict{Point{S}, Int}) where S <: FieldElem

Given a curve `C` which is assumed to be smooth and irreducible, return the divisor on the curve `C` defined by `D`.
"""
struct AffineCurveDivisor{S <: FieldElem} <: CurveDivisor
    C::AffinePlaneCurve{S}
    divisor::Dict{Point{S}, Int}
    degree::Int
    function AffineCurveDivisor(C::AffinePlaneCurve{S}, D::Dict{Point{S}, Int}) where S <: FieldElem
        for (P, m) in D
            P in C || error("The point ", P.coord, " is not on the curve")
            if m == 0
                delete!(D, P)
            end
        end
        new{S}(C, D, sum(values(D)))
    end
end

function AffineCurveDivisor(C::AffinePlaneCurve{S}, P::Point{S}, m::Int=1) where S <: FieldElem
    return AffineCurveDivisor(C, Dict(P => m))
end

################################################################################

@doc Markdown.doc"""
    ProjCurveDivisor(C::ProjPlaneCurve{S}, D::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}) where S <: FieldElem

Given a curve `C` which is assumed to be smooth and irreducible, return the divisor on the curve `C` defined by `D`.
"""
struct ProjCurveDivisor{S <: FieldElem} <: CurveDivisor
    C::ProjPlaneCurve{S}
    divisor::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}
    degree::Int
    function ProjCurveDivisor(C::ProjPlaneCurve{S}, D::Dict{Oscar.Geometry.ProjSpcElem{S}, Int}) where S <: FieldElem
        for (P, m) in D
            P in C || error("The point ", P, " is not on the curve")
            if m == 0
                delete!(D, P)
            end
        end
        new{S}(C, D, sum(values(D)))
    end
end

function ProjCurveDivisor(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}, m::Int=1) where S <: FieldElem
    return ProjCurveDivisor(C, Dict(P => m))
end

################################################################################
# To get a nice output for the divisors. The needs_parentheses is needed since
# expressify wants to know if parentheses are needed.

function AbstractAlgebra.needs_parentheses(A::Array{T, 1}) where T <: RingElem
    return false
end

function AbstractAlgebra.expressify(@nospecialize(a::AffineCurveDivisor); context = nothing)
   prod = Expr(:call, :+)
   for (P, m) in a.divisor
      ep = AbstractAlgebra.expressify(P.coord, context = context)
      push!(prod.args, Expr(:call, :*, m, ep))
   end
   return prod
end

function AbstractAlgebra.needs_parentheses(A::Oscar.Geometry.ProjSpcElem{T}) where T <: RingElem
    return false
end

function AbstractAlgebra.expressify(@nospecialize(a::ProjCurveDivisor); context = nothing)
   prod = Expr(:call, :+)
   for (P, m) in a.divisor
      ep = AbstractAlgebra.expressify(P, context = context)
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

function curve(D::CurveDivisor)
    return D.C
end

function Oscar.degree(D::CurveDivisor)
    return D.degree
end

################################################################################

function _check_same_curve(D::CurveDivisor, E::CurveDivisor)
    D.C == E.C || error("The divisors are not on the same curve")
end

################################################################################
################################################################################
# Sum of two divisors

function Base.:+(D::AffineCurveDivisor, E::AffineCurveDivisor)
    _check_same_curve(D, E)
    return AffineCurveDivisor(D.C, merge(+, D.divisor, E.divisor))
end

function Base.:+(D::ProjCurveDivisor, E::ProjCurveDivisor)
    _check_same_curve(D, E)
    return ProjCurveDivisor(D.C, merge(+, D.divisor, E.divisor))
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

function Base.:-(D::AffineCurveDivisor, E::AffineCurveDivisor)
    return AffineCurveDivisor(D.C, _help_minus(D, E))
end

function Base.:-(D::ProjCurveDivisor, E::ProjCurveDivisor)
    return ProjCurveDivisor(D.C, _help_minus(D, E))
end

################################################################################
################################################################################

function Base.:*(k::Int, D::AffineCurveDivisor{S}) where S <: FieldElem
    if k == 0
        return AffineCurveDivisor(D.C, Dict{Point{S}, Int}())
    elseif k == 1
        return D
    else
        return AffineCurveDivisor(D.C, Dict((P, k*m) for (P, m) in D.divisor))
    end
end

function Base.:*(k::Int, D::ProjCurveDivisor{S}) where S <: FieldElem
    if k == 0
        return ProjCurveDivisor(D.C, Dict{Oscar.Geometry.ProjSpcElem{S}, Int}())
    elseif k == 1
        return D
    else
        return ProjCurveDivisor(D.C, Dict((P, k*m) for (P, m) in D.divisor))
    end
end

################################################################################

function ==(D::CurveDivisor, E::CurveDivisor)
    return D.C == E.C && D.divisor == E.divisor
end

################################################################################
################################################################################

@doc Markdown.doc"""
    iseffective(D::CurveDivisor)

Return `true` if `D` is an effective divisor, `false` otherwise.
"""
function iseffective(D::CurveDivisor)
    return all(v >= 0 for v in values(D.divisor))
end

################################################################################
################################################################################
# TODO: take phi in the function field of C when fixed

@doc Markdown.doc"""
    multiplicity(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Point{S}) where {S <: FieldElem, T <: MPolyElem{S}}

Return the multiplicity of the rational function `phi` on the curve `C` at the point `P`.
"""
function multiplicity(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Point{S}) where {S <: FieldElem, T <: MPolyElem{S}}
    f = divrem(phi.num, C.eq)
    g = divrem(phi.den, C.eq)
    !iszero(g[2]) || error("This is not a rational function on `C`")
    return multiplicity(C, f[2], P) - multiplicity(C, g[2], P)
end

###############################################################################

@doc Markdown.doc"""
    multiplicity(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}, P::Point{S}) where S <: FieldElem

Return the multiplicity of the polynomial `F` on the curve `C` at the point `P`.
"""
function multiplicity(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}, P::Point{S}) where S <: FieldElem
    f = divrem(F, C.eq)
    if iszero(f[2])
        return Inf
    elseif isconstant(f[2])
        return 0
    else
        return intersection_multiplicity(C, AffinePlaneCurve(f[2]), P)
    end
end

################################################################################

@doc Markdown.doc"""
    multiplicity(C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the multiplicity of the polynomial `F` on the curve `C` at the point `P`.
"""
function multiplicity(C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
    f = divrem(F.f, defining_equation(C))
    if iszero(f[2])
        return Inf
    elseif isconstant(f[2])
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
    if !isconstant(f[2])
        CC = AffinePlaneCurve(f[2])
        L = curve_intersect(C, CC)
        for P in L[2]
            E[P] = intersection_multiplicity(C, CC, P)
        end
    end
    return AffineCurveDivisor(C, E)
end

################################################################################
# TODO: take phi in the function field of C when fixed

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
    if !isconstant(f[2])
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
    PP = projective_space(R.R.base_ring, 2)
    return divisor(PP[1], C, F)
end

################################################################################
# TODO divisor with fraction of decorated polynomials when available
