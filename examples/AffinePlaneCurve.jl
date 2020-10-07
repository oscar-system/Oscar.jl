module AffinePlaneCurveModule

using Oscar, Markdown

export AffinePlaneCurve, degree, jacobi_ideal, tangent, issmooth_point, union,
       Point

################################################################################
#
# We follow the lecture notes by Andreas Gathmann: 
# https://www.mathematik.uni-kl.de/~gathmann/en/curves.php
#
################################################################################
mutable struct AffinePlaneCurve{S}
  eq::Oscar.MPolyElem{S}                  # Equation of the curve (polynomial in two variables)
  degree::Int64                           # degree of the equation of the curve
  dimension::Int64                        # dimension of the curve (as a variety)
  function AffinePlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    r = new{S}()
    r.eq = eq
    if nvars(parent(r.eq)) != 2           
       error("The defining equation must belong to a ring with two variables")
    else
       r.degree = -1                      # -1 when it is not computed yet
       r.dimension = 1                    # since C is a plane curve, the dimension is always 1
       return r
    end
  end
end

function Base.show(io::IO, C::AffinePlaneCurve)
  println(io, "Affine plane curve defined by ", C.eq)
end

################################################################################
# Point (first attempt).
# The point is described by its coordinates.
# ambient_dim gives the dimension of the space to which belongs the point.

mutable struct Point{S}
  coord::Array{S, 1}
  ambient_dim::Int64
  function Point(coordinates::Array{S,1}) where {S <: FieldElem}
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
# Compute the degree of the equation of the curve if not already known, and show it.

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
# To check if a point is smooth: return true if the point is a smooth point on the curve C, and false  if it is singular, and an error if it is not on the curve.

@doc Markdown.doc"""
issmooth_point(C::AffinePlaneCurve{S}, P::Point{S}) where S<: FieldElem
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
     a = evaluate(L[1],P.coord)
     b = evaluate(L[2],P.coord)
     if a == 0 && b == 0 
        return false
     else
        return true
     end
  end
end

################################################################################
# If P is a smooth point of an affine plane curve C, compute the tangent of C at the point P.

@doc Markdown.doc"""
tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S<: FieldElem
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
     a = evaluate(L[1],P.coord)
     b = evaluate(L[2],P.coord)
     if a == 0 && b == 0
        error("This is a singular point of the curve")
     else
        V = gens(parent(C.eq))
        T = a*( V[1] - P.coord[1] ) + b*( V[2] - P.coord[2] )
        return AffinePlaneCurve(T)
     end
  end
end

################################################################################
# Union of two affine plane curves (with multiplicity)

@doc Markdown.doc"""
union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S<: FieldElem
> Given two affine plane curves C and D, returns the union of the two curves (with multiplicity).
"""
function Base.union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  F = C.eq*D.eq
  return AffinePlaneCurve(F)
end

################################################################################
end
using .AffinePlaneCurveModule
