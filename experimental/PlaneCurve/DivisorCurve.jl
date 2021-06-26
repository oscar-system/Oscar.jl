export AffineCurveDivisor, ProjCurveDivisor, iseffective, multiplicity,
       divisor, curve, degree, divisor_ideals, global_sections,
       dimension_global_sections, islinearly_equivalent, isprincipal,
       curve_zero_divisor, principal_divisor

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
# To get a nice output for the divisors. The needs_parentheses is needed since
# expressify wants to know if parentheses are needed.

function AbstractAlgebra.needs_parentheses(A::Vector{T}) where T <: RingElem
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

function curve_zero_divisor(C::ProjPlaneCurve{S}) where S <: FieldElem
    return ProjCurveDivisor{S}(C, Dict{Oscar.Geometry.ProjSpcElem{S}, Int}())
end

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
    iseffective(D::CurveDivisor)

Return `true` if `D` is an effective divisor, `false` otherwise.
"""
function iseffective(D::CurveDivisor)
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

function multiplicity(C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Oscar.Geometry.ProjSpcElem{S})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
    g = divrem(phi.den.f, C.eq.f)
    !iszero(g[2]) || error("This is not a rational function on the curve")
    f = divrem(phi.num.f, C.eq.f)
    !iszero(f[2]) || error("The numerator is zero on the curve")
    R = parent(C.eq)
    return multiplicity(C, R(f[2]), P) - multiplicity(C, R(g[2]), P)
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
    PP = projective_space(R.R.base_ring, 2)
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

# colon ideal (I : J) modulo Q
function _qquotient(I::Oscar.MPolyIdeal, J::Oscar.MPolyIdeal, Q::Oscar.MPolyIdeal)
   R = base_ring(I)
   R == base_ring(J) || error("Cannot be computed")
   R == base_ring(Q) || error("Cannot be computed")
   II = I + Q
   JJ = J + Q
   Q1 = quotient(II, JJ)
   B = groebner_basis(Q, ordering = :lex, complete_reduction = true)
   r = divrem(gens(Q1), B)
   arr = [r[i][2] for i in 1:length(r)]
   return ideal(R, arr)
end

function jet(I::Oscar.MPolyIdeal, d::Int)
   singular_assure(I)
   R = base_ring(I)
   Is = Singular.jet(I.gens.S, d)
   return ideal(R, Is)
end

function _remove_zeros(I::Oscar.MPolyIdeal)
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
    return ideal(R, Imin)
end

# Remove components which are not codim 1
function _purify1(I::Oscar.MPolyIdeal, Q::Oscar.MPolyIdeal)
   R = base_ring(I)
   Id = _remove_zeros(I)
   gens(Id)[1] != R(0) || error("ideal assumed to be non-zero")
   f = ideal(R, [gens(Id)[1]])
   return _minimal_generating_set(_qquotient(f, _qquotient(f, Id, Q), Q))
end

function _basis(I::Oscar.MPolyIdeal, d::Int)
   R = base_ring(I)
   m = ideal(R, gens(R))
   return _minimal_generating_set(_remove_zeros(jet(intersect(I, m^d),d)))
end

function _global_sections_helper(I::Oscar.MPolyIdeal, J::Oscar.MPolyIdeal, Q::Oscar.MPolyIdeal)
   R = base_ring(I)
   R == base_ring(J) || error("base rings do not match")
   R == base_ring(Q) || error("base rings do not match")
   Ip = _purify1(I, Q)
   Jp = _purify1(J, Q)
   Is = _remove_zeros(Ip)
   f = gens(Is)[1]
   fJ = f*Jp
   Q1 = _qquotient(fJ, Ip, Q)
   P = _purify1(Q1, Q)
   B = _basis(P, total_degree(f))
   BQ = groebner_basis(Q, ordering = :lex, complete_reduction = true)
   r = divrem(gens(B), BQ)
   arr = [r[i][2] for i in 1:length(r)]
   LD = ideal(R, arr)
   return [_remove_zeros(LD), f]
end

################################################################################

function _global_sections_ideals(D::ProjCurveDivisor)
    F = defining_equation(D.C)
    R = parent(F)
    Q = ideal(R, [F])
    E = divisor_ideals(D)
    _global_sections_helper(E[1], E[2], Q)
end

################################################################################

@doc Markdown.doc"""
    global_sections(D::ProjCurveDivisor)
Return a set of generators of the global sections of the sheaf associated to the divisor `D` of a smooth and irreducible projective curve.
"""
function global_sections(D::ProjCurveDivisor)
    L = _global_sections_ideals(D)
    N = gens(L[1])
    return [g//L[2] for g in N]
end

################################################################################

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
    J = ideal(R, [g])
    Q = ideal(R, [F])
    return _purify1(_qquotient(_qquotient(f*E[1], J, Q), E[2], Q), Q)
end

################################################################################

function _linearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
    F = D - E
    R = parent(defining_equation(D.C))
    T = parent(D.C.eq)
    L = _global_sections_ideals(F)
    if length(gens(L[1])) !=1
        return T(0)//T(1)
    elseif iszero(gens(L[1])[1])
        return T(0)//T(1)
    else
        V = _section_ideal(gens(L[1])[1], L[2], F)
        if V == ideal(R, [R(1)])
            return T(gens(L[1])[1])//T(L[2])
        else
            return T(0)//T(1)
        end
    end
end

################################################################################
@doc Markdown.doc"""
    islinearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
Return `true` if the divisors `D` and `E` are linearly equivalent, and `false`
otherwise"""
function islinearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
    return !iszero(_linearly_equivalent(D, E))
end

################################################################################
@doc Markdown.doc"""
    isprincipal(D::ProjCurveDivisor{S}) where S <: FieldElem
Return `true` if the divisor `D` is principal, and `false` otherwise"""
function Oscar.isprincipal(D::ProjCurveDivisor{S}) where S <: FieldElem
    if !isdefined(D, :prin_div)
        D.prin_div = _linearly_equivalent(D, curve_zero_divisor(D.C))
    end
    return !iszero(D.prin_div)
end

################################################################################
@doc Markdown.doc"""
    principal_divisor(D::ProjCurveDivisor{S}) where S <: FieldElem
If the divisor `D` is principal, return a rational function `phi` such that `D`
is linearly equivalent to the divisor defined by `phi`."""
function principal_divisor(D::ProjCurveDivisor{S}) where S <: FieldElem
    isprincipal(D) || error("the divisor is not principal")
    return D.prin_div
end
################################################################################
