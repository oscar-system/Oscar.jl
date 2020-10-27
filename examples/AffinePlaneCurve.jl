module AffinePlaneCurveModule

using Oscar, Markdown
include("Variety.jl")
using .VarietyModule
export AffinePlaneCurve, degree, jacobi_ideal, tangent, issmooth_point, union,
       Point, curve_components, isirreducible, isreduced, common_components,
       curve_intersect, intersect, curve_singular_locus

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
       r.components = Dict()
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
> Given an affine plane curve C, returns the degree of the defining polynomial.
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
> Given an affine plane curve C, returns the Jacobian ideal of the defining polynomial.
"""
function Oscar.jacobi_ideal(C::AffinePlaneCurve)
 return jacobi_ideal(C.eq)
end

################################################################################
# To check if a point is smooth: return true if the point is a smooth point on 
# the curve C, and false  if it is singular, and an error if it is not on the 
# curve.

@doc Markdown.doc"""
    issmooth_point(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Given an affine plane curve C and a point P, returns an error if the point is not in the zero locus of the defining equation, false if it is a singular point of C, and true if it is a smooth point of the curve.
"""
function issmooth_point(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
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
> Given an affine plane curve C and a point P, if P is a smooth point of C, returns the affine plane curve defined as the tangent of C at P.
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
> Given two affine plane curves C and D, returns the union of the two curves (with multiplicity).
"""
function Base.union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  F = C.eq*D.eq
  return AffinePlaneCurve(F)
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
    curve_components(C::AffinePlaneCurve)
> Given an affine plane curve C, returns a dictionary containing its irreducible components (as affine plane curves) and their multiplicity.
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
> Given an affine plane curve C, returns true if the curve is irreducible, and false otherwise.
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

@doc Markdown.doc"""
    isreduced(C::AffinePlaneCurve)
> Given an affine plane curve C, returns true if the curve is reduced, and false otherwise.
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
    common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem(C::AffinePlaneCurve)
> Given two affine plane curves C and D, returns the affine plane curve consisting of the common components, or an error if they do not have a common component.
"""
function common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  if G == 1
     error("The curves do not have a common component")
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
> Given two affine plane curves C and D defined by F and G, returns a list whose first element is the affine plane curve defined by the gcd of F and G, the second element is the list of the remaining intersection points when the common components are removed from C and D.
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
> Given two affine plane curves C and D, returns the variety defined by the intersection of the two curves
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
> Given an affine plane curve C, returns a list whose first element is the affine plane curve consisting of the singular components of C (if any), and the second element is the list of the isolated singular points (which may be contained in the singular component). The singular component might not contain any point over the considered field.
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
           return L
        end
     end
  else
     if isconstant(FY) == true
        if FY != 0
           return [CC, Pts]
        else
           CX = AffinePlaneCurve(FX)
           L = curve_intersect(C, CX)
           return L
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
           CC = M[1]
           PP = [Pts; M[2]]
           return [CC, PP]
        end
     end
  end
end

################################################################################

end
using .AffinePlaneCurveModule
