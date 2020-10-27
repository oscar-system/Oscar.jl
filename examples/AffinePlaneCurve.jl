module AffinePlaneCurveModule

using Oscar, Markdown
include("Variety.jl")
using .VarietyModule

import Base.:(==)

export AffinePlaneCurve, degree, jacobi_ideal, tangent, issmooth, union,
       Point, curve_components, isirreducible, isreduced, common_component,
       curve_intersect, intersect, curve_singular_locus, ideal_point,
       multiplicity, singular_locus, intersection_multiplicity,
       aretransverse, tangent_lines, issmooth_curve, reduction, hash

################################################################################
#
# We follow the lecture notes by Andreas Gathmann: 
# https://www.mathematik.uni-kl.de/~gathmann/en/curves.php
#
################################################################################
#

mutable struct AffinePlaneCurve{S}
  eq::Oscar.MPolyElem{S}                  # Equation of the curve (polynomial in two variables)
  degree::Int64                           # degree of the equation of the curve
  dimension::Int64                        # dimension of the curve (as a variety)
  components::Dict
  function AffinePlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    r = new{S}()
    r.eq = eq
    if nvars(parent(r.eq)) != 2           
       error("The defining equation must belong to a ring with two variables")
    elseif isconstant(r.eq) == true
       error("The defining equation must be non constant")
    else
       r.degree = -1                      # -1 when it is not computed yet
       r.dimension = 1                    # since C is a plane curve, the dimension is always 1
       r.components = Dict{AffinePlaneCurve{S}, Int64}()
       return r
    end
  end
end

function Base.show(io::IO, C::AffinePlaneCurve)
  if !get(io, :compact, false)
     println(io, "Affine plane curve defined by ", C.eq)
  else
     println(io, C.eq)
  end
end

################################################################################
# Point (first attempt).
# The point is described by its coordinates.
# ambient_dim gives the dimension of the space to which belongs the point.

mutable struct Point{S}
  coord::Array{S, 1}
  ambient_dim::Int64
  function Point(coordinates::Array{S, 1}) where {S <: FieldElem}
    r = new{S}()
    r.coord = coordinates
    r.ambient_dim = length(coordinates)
    return r
  end
end

function Base.show(io::IO, P::Point)
  println(io, "Point with coordinates ", P.coord)
end

################################################################################
# Compute the degree of the equation of the curve if not already known, 
# and show it.

@doc Markdown.doc"""
degree(C::AffinePlaneCurve)
> Returns the degree of the defining polynomial of C.
"""
function Oscar.degree(C::AffinePlaneCurve)
  if C.degree >= 0
     return C.degree
  else
     C.degree = total_degree(C.eq)
     return(C.degree)
  end
end

################################################################################
# Compute the Jacobian Ideal of C

@doc Markdown.doc"""
jacobi_ideal(C::AffinePlaneCurve)
> Returns the Jacobian ideal of the defining polynomial of C.
"""
function Oscar.jacobi_ideal(C::AffinePlaneCurve)
 return jacobi_ideal(C.eq)
end

################################################################################
# To check if a point is smooth: return true if the point is a smooth point on 
# the curve C, and false  if it is singular, and an error if it is not on the 
# curve.

@doc Markdown.doc"""
issmooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Returns an error if P is not a point of C, false if P is a singular point of C, and true if P is a smooth point of C.
"""
function Oscar.issmooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  if P.ambient_dim != 2
     error("The point needs to be in a two dimensional space")
  elseif evaluate(C.eq, P.coord) != 0
     error("The point is not on the curve defined by ", C.eq)
  else
     J = jacobi_ideal(C)
     L = gens(J)
     a = evaluate(L[1], P.coord)
     b = evaluate(L[2], P.coord)
     if a == 0 && b == 0 
        return false
     else
        return true
     end
  end
end

################################################################################
# If P is a smooth point of an affine plane curve C, compute the tangent of C 
# at the point P.

@doc Markdown.doc"""
tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Returns the tangent of C at P when P is a smooth point of C, and an error otherwise.
"""
function tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  if P.ambient_dim != 2
     error("The point needs to be in a two dimensional space")
  elseif evaluate(C.eq, P.coord) != 0
     error("The point is not on the curve")
  else
     J = jacobi_ideal(C)
     L = gens(J)
     a = evaluate(L[1], P.coord)
     b = evaluate(L[2], P.coord)
     if a == 0 && b == 0
        error("This is a singular point of the curve")
     else
        V = gens(parent(C.eq))
        T = a*(V[1] - P.coord[1]) + b*(V[2] - P.coord[2])
        return AffinePlaneCurve(T)
     end
  end
end

################################################################################
# Union of two affine plane curves (with multiplicity)

@doc Markdown.doc"""
union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
> Returns the union of C and D (with multiplicity).
"""
function Base.union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  F = C.eq*D.eq
  return AffinePlaneCurve(F)
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
curve_components(C::AffinePlaneCurve)
> Returns a dictionary containing the irreducible components of C and their multiplicity.
"""
function curve_components(C::AffinePlaneCurve)
  if C.components != Dict()
     return C.components
  else
     D = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
     return C.components
  end
end

################################################################################
# Check irreducibility.
# TODO: change for a direct irreducibility check when available.

@doc Markdown.doc"""
isirreducible(C::AffinePlaneCurve)
> Returns true if C is irreducible, and false otherwise.
"""
function Oscar.isirreducible(C::AffinePlaneCurve)
  if C.components == Dict()
     D = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  if length(C.components) != 1
     return false
  elseif collect(values(C.components)) != [1]
     return false
  else
     return true
  end
end

################################################################################
# Check reducedness by computing a factorization
# TODO: change for a direct squarefree check when available.

@doc Markdown.doc"""
isreduced(C::AffinePlaneCurve)
> Returns true if C is reduced, and false otherwise.
"""
function Oscar.isreduced(C::AffinePlaneCurve)
  if C.components == Dict()
     D = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  if [C.components[x] for x in keys(C.components)] == [1 for x in keys(C.components)]
     return true
  else
     return false
  end
end

################################################################################
# gives the common components of two affine plane curves

@doc Markdown.doc"""
common_component(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem(C::AffinePlaneCurve)
> Returns the affine plane curve consisting of the common component of C and D, or an empty array if they do not have a common component.
"""
function common_component(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  if G == 1
     return(Array{AffinePlaneCurve, 1}())
  else
     return AffinePlaneCurve(G)
  end
end

################################################################################
# Compute the intersection of two curves. The output is a list: the first
# element the affine plane curve defined by the common components of the two
# curves (or is empty if no common component), the second element is the list
# of intersection points. Some might be contained in the common component too.

@doc Markdown.doc"""
curve_intersect(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
> Returns a list whose first element is the common component of C and D, and the second element is the list of the remaining intersection points.
"""
function curve_intersect(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  R = parent(C.eq)
  if G != 1
     # We divide by the gcd to get curves without common components.
     F = div(C.eq, G)
     H = div(D.eq, G)
  else
     F = C.eq
     H = D.eq
  end
  I = ideal(R, [F, H])
  # We eliminate the first variable in the ideal.
  B = eliminate(I, [gens(R)[1]])
  # We factorize the polynomial we obtain (which is a polynomial where only
  # the second variable appears).
  Z = factor(gens(B)[1])
  Y = []
  #L = []
  L = Array{S, 1}[]
  # For the linear factors of g, we add the constant coefficient to the list Y.
  # The constant coefficient is stored as a mpoly to be used later in evaluate.
  for g in keys(Z.fac)
     if total_degree(g) == 1
        f = g//lc(g)
        push!(Y, -f + gens(R)[2])
     end
  end
  if Y != []
     # For each y, we compute the possible values of x by replacing the second
     # variable by the value y, and factorizing the resulting polynomial.
     for y in Y
        FF = evaluate(F, [gens(R)[1], y])
        HH = evaluate(H, [gens(R)[1], y])
        GG = gcd(FF, HH)
        ZZ = factor(GG)
        for g in keys(ZZ.fac)
           if total_degree(g) == 1
              f = g//lc(g)
              # We use lc to convert the mpoly into element of S.
              push!(L, [lc(-f + gens(R)[1]), lc(y)])
           end
        end
     end
  end
  # L contains the intersection points as arrays (or can be empty).
  Pts = Array{Point, 1}()
  CC = Array{AffinePlaneCurve, 1}()
  if L != []
     for p in L
        push!(Pts, Point(p))
     end
  end
  if G != 1
     push!(CC, AffinePlaneCurve(G))
  end
  if L == [] && G == 1
     println("The curves do not intersect")
     return [CC, Pts]
  else
     return [CC, Pts]
  end
end

################################################################################
# Intersection in the sense of varieties.
# Might change depending on the futur changes on VarietyModule.

@doc Markdown.doc"""
intersect( C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
> Returns the variety defined by the intersection of C and D.
"""
function Oscar.intersect(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  I = ideal(parent(C.eq), [C.eq, D.eq])
  V = Variety(I)
  return V
end

################################################################################
# Compute the singular locus of an affine plane curve by intersecting the
# derivatives and then checking if the intersection belongs also to the
# curve. The points might also be contained in the components.

@doc Markdown.doc"""
curve_singular_locus( C::AffinePlaneCurve)
> Returns a list whose first element is the affine plane curve consisting of the singular components of C (if any), and the second element is the list of the isolated singular points of C (which may be contained in the singular component). The singular component might not contain any point over the considered field.
"""
function curve_singular_locus(C::AffinePlaneCurve)
  F = C.eq
  J = jacobi_ideal(F)
  FX = gens(J)[1]
  FY = gens(J)[2]
  Pts = Array{Point, 1}()
  CC = Array{AffinePlaneCurve, 1}()
  # Case where some derivatives are constant. If one of them is constant
  # non-zero, the singular locus is empty. If both derivatives are equal to 0,
  # the curve is singular at all its points. If one of them is 0 and the other
  # one non constant, the intersection of C and the non constant derivative
  # gives the singular locus.
  if isconstant(FX) == true
     if FX != 0
        return [CC, Pts]
     else
        if isconstant(FY) == true
           if FY == 0
              push!(CC, C)
              return [CC, Pts]
           else
              return [CC, Pts]
           end
        else
           CY = AffinePlaneCurve(FY)
           L = curve_intersect(C, CY)
           M = L[1]
           rCC = reduction(M[1])
           push!(CC, AffinePlaneCurve(rCC.eq*M[1].eq))
           return [CC, L[2]]
        end
     end
  else
     if isconstant(FY) == true
        if FY != 0
           return [CC, Pts]
        else
           CX = AffinePlaneCurve(FX)
           L = curve_intersect(C, CX)
           M = L[1]
           rCC = reduction(M[1])
           push!(CC, AffinePlaneCurve(rCC.eq*M[1].eq))
           return [CC, L[2]]
        end
     else
        # General case: both derivatives are non constant. We intersect the
        # curves defined by the derivatives, and then check the intersection
        # with the curve.
        CX = AffinePlaneCurve(FX)
        CY = AffinePlaneCurve(FY)
        S = curve_intersect(CX, CY)
        for p in S[2]
           if evaluate(F, p.coord) == 0
              push!(Pts, p)
           end
        end
        if S[1] == []
           return [CC, Pts]
        else
           M = curve_intersect(C, S[1][1])
           MM = M[1]
           PP = [Pts; M[2]]
           rCC = reduction(MM[1])
           push!(CC, AffinePlaneCurve(rCC.eq*MM[1].eq))
           return [CC, PP]
        end
     end
  end
end

################################################################################
# Compute squarefree equation.
# TODO: change for a direct squarefree computation when available.

@doc Markdown.doc"""
reduction(C::AffinePlaneCurve)
> Returns the affine plane curve defined by the squarefree part of the equation of C.
"""
function reduction(C::AffinePlaneCurve)
  if C.components == Dict()
     f = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => f.fac[x] for x in keys(f.fac))
  end
  F = 1
  for D in keys(C.components)
     F = F*D.eq
  end
  rC = AffinePlaneCurve(F)
  rC.components = Dict(AffinePlaneCurve(D.eq) => 1 for D in keys(C.components))
  return(rC)
end


################################################################################
# Associate a maximal ideal to a point in a given ring (not specific to curves)

@doc Markdown.doc"""
ideal_point(P::Point{S}, R::MPolyRing{S}) where S <: FieldElem
> Returns the maximal ideal associated to the point P in the ring R.
"""
function ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem
  V = gens(R)
  return ideal(R, [V[i] - P.coord[i] for i in 1:length(V)])
end

################################################################################
# Helping function:
# change of variable to send a point at the origin.

function curve_map_point_origin(C::AffinePlaneCurve, P::Point)
  F = C.eq
  R = parent(C.eq)
  V = gens(R)
  G = evaluate(F, [V[1] + P.coord[1], V[2] + P.coord[2]])
  return(AffinePlaneCurve(G))
end

################################################################################
# Helping function:
# used to order elements to find the multiplicity.

function _sort_helper_multiplicity(a::GrpAbFinGenElem)
  return a.coeff[1, 1]
end

################################################################################
# compute the multiplicity of the affine plane curve C at the point P

@doc Markdown.doc"""
multiplicity(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Returns the multiplicity of C at P.
"""
function multiplicity(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  if P.ambient_dim != 2
     error("The point needs to be in a two dimensional space")
  end
  D = curve_map_point_origin(C, P)
  G = D.eq
  R = parent(G)
  A = grade(R)
  HC = homogenous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  return(M[1].coeff[1, 1])
end

################################################################################
# compute the set of tangent lines of the affine plane curve C at the point P
# (linear factors of the homogeneous part of lower degree of the equation).

@doc Markdown.doc"""
tangent_lines(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Returns the tangent lines at P to C with their multiplicity.
"""
function tangent_lines(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  if P.ambient_dim != 2
     error("The point needs to be in a two dimensional space")
  end
  D = curve_map_point_origin(C, P)
  G = D.eq
  R = parent(G)
  V = gens(R)
  A = grade(R)
  HC = homogenous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  Gm = HC[M[1]]
  Z = factor(Gm)
  D = Dict{AffinePlaneCurve{S}, Int64}()
  X = V[1] - P.coord[1]
  Y = V[2] - P.coord[2]
  for p in keys(Z.fac)
     if total_degree(p.f) == 1
        push!(D, AffinePlaneCurve(evaluate(p.f, [X, Y])) => Z.fac[p])
     end
  end
  return(D)
end

################################################################################
# check for equality of curves (the equations are equal up to multiplication by
# a non zero constant).

function ==(C::AffinePlaneCurve, D::AffinePlaneCurve)
  F = C.eq
  G = D.eq
  d = div(F, G)
  if isconstant(d) == false
     return false
  else
     if F - d*G == 0
        return true
     else
        return false
     end
  end
end

################################################################################
# Singular locus in the sense of varieties.
# Might change depending on the futur changes on VarietyModule.

@doc Markdown.doc"""
singular_locus(C::AffinePlaneCurve)
> Returns the singular locus of C as a variety.
"""
function singular_locus(C::AffinePlaneCurve)
  J = jacobi_ideal(C)
  R = parent(C.eq)
  S = ideal(R, [gens(J); C.eq])
  return Variety(S)
end

################################################################################
# Compute the intersection multiplicity of the two affine plane curves at the
# given point.

@doc Markdown.doc"""
intersection_multiplicity(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Returns the intersection multiplicity of C and D at P.
"""
function intersection_multiplicity(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  if P.ambient_dim != 2
     error("The point needs to be in a two dimensional space")
  end
  R = parent(C.eq)
  m = ideal_point(R, P)
  r = Localization(R, m)
  I = ideal(r, [C.eq, D.eq])
  groebner_basis(I)
  return(Singular.vdim(I.gb.S))
end

################################################################################
# Check if the two curves intersect transversally at the given point.

@doc Markdown.doc"""
aretransverse(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S<:FieldElem
> Returns true if C and D intersect transversally at P and false otherwise.
"""
function aretransverse(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S<:FieldElem
  if P.ambient_dim != 2
     error("The point needs to be in a two dimensional space")
  end
  if issmooth(C, P) == false || issmooth(D, P) == false
     return false
  else
     if intersection_multiplicity(C, D, P) == 1
        return true
     else
        return false
     end
  end
end

################################################################################
# Check if a reduced curve is smooth.

@doc Markdown.doc"""
issmooth_curve(C::AffinePlaneCurve)
> Given a reduced affine plane curve C, returns true if C has no singular point, and false otherwise.
"""
function issmooth_curve(C::AffinePlaneCurve)
  S = curve_singular_locus(C)
  if S[1] == []
     if S[2] == []
        return true
     else
        return false
     end
  else
     error("The curve is not reduced.")
  end
end

################################################################################
# hash functions

function Base.hash(C::AffinePlaneCurve, h::UInt)
  F = 1//lc(C.eq)*C.eq
  return hash(F, h)
end

################################################################################

end
using .AffinePlaneCurveModule
