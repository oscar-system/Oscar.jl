module AffinePlaneCurveModule

using Oscar, Markdown

export AffinePlaneCurve, AffinePlaneCurve_degree, AffinePlaneCurve_Jacobian, AffinePlaneCurve_tangent, AffinePlaneCurve_issmoothpoint, AffinePlaneCurve_union

####################################################
#
#We follow the lecture notes by Andreas Gathmann: https://www.mathematik.uni-kl.de/~gathmann/en/curves.php
#
####################################################
mutable struct AffinePlaneCurve{S}
  eq::Oscar.MPolyElem                      #Equation of the curve (polynomial in two variables)
  degree::Union{String, Int64}             #degree of the equation of the curve
  dimension::Int64                         #dimension of the curve (as a variety)
  function AffinePlaneCurve(eq::Oscar.MPolyElem{T}) where {T <: FieldElem}
    r = new{T}()
    r.eq = eq
    if nvars(parent(r.eq))!=2              #If the ambiant ring does not have two variables, show an error.
       error("The defining equation must belong to a ring with two variables")
    else
       r.degree = "not computed yet"
       r.dimension = 1                    #since C is a plane curve, the dimension is always 1
       return r
    end
  end
end

function Base.show(io::IO, C::AffinePlaneCurve)
  println(io, "Affine plane curve defined by ", C.eq)
end

######################################################
#Compute the degree of the equation of the curve if not already known, and show it.


@doc Markdown.doc"""
   AffinePlaneCurve_degree(C::AffinePlaneCurve)
   > Given an affine plane curve C, returns the degree of the defining polynomial.
"""
function AffinePlaneCurve_degree(C::AffinePlaneCurve)
  if typeof(C.degree)==Int64
     println(C.degree)
  else
  C.degree=total_degree(C.eq)
  println(C.degree)
  end
end

######################################################
#Compute the Jacobian Ideal of C

@doc Markdown.doc"""
   AffinePlaneCurve_Jacobian(C::AffinePlaneCurve)
   > Given an affine plane curve C, returns the Jacobian ideal of the defining polynomial.
"""
function AffinePlaneCurve_Jacobian(C::AffinePlaneCurve)
 return jacobi_ideal(C.eq)
end

######################################################
#To check if a point is smooth: return true if the point is a smooth point on the curve C, and false otherwise.

@doc Markdown.doc"""
   AffinePlaneCurve_Jacobian(C::AffinePlaneCurve, P::Array{S,1}) where S<: FieldElem
   > Given an affine plane curve C and a point P, returns false if the point is not in the zero locus of the defining equation or is a singular point of C, and true if it is a smooth point of the curve.
"""
function AffinePlaneCurve_issmoothpoint(C::AffinePlaneCurve{S}, P::Array{S,1}) where S<: FieldElem
  if evaluate(C.eq, P)!=0
     return false                         #return false if the point is not on the curve
  else
     J=AffinePlaneCurve_Jacobian(C)
     L=gens(J)
     a=evaluate(L[1],P)
     b=evaluate(L[2],P)
     if a==0 && b==0 
        return false                      #return false if the point is not a smooth point of the curve
     else
        return true
     end
  end
end

######################################################
#Compute the tangent of C at the point (x0,y0), if this point is a smooth point of C
#For a point over QQ, it has to be written as P = [QQ(1), QQ(-1//2)]

@doc Markdown.doc"""
   AffinePlaneCurve_tangent(C::AffinePlaneCurve, P::Array{S,1}) where S<: FieldElem
   > Given an affine plane curve C and a point P, if P is a smooth point of C, returns the affine plane curve defined as the tangent of C at P.
"""
function AffinePlaneCurve_tangent(C::AffinePlaneCurve{S}, P::Array{S, 1}) where S <: FieldElem
  if evaluate(C.eq, P)!=0
     println("The point is not on the curve")
  else
     J=AffinePlaneCurve_Jacobian(C)
     L=gens(J)
     a=evaluate(L[1],P)
     b=evaluate(L[2],P)
     if a==0 && b==0
        println("This is not a smooth point of the curve")
     else
        V=gens(parent(C.eq))
        T=a*(V[1]-P[1])+b*(V[2]-P[2])
        return AffinePlaneCurve(T)
     end
  end
end

####################################################
#Union of two affine plane curves (with multiplicities)

@doc Markdown.doc"""
   AffinePlaneCurve_union(C::AffinePlaneCurve, D::AffinePlaneCurve{S}) where S<: FieldElem
   > Given two affine plane curves C and D, returns the union of the two curves (with multiplicity).
"""
function AffinePlaneCurve_union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  F=C.eq*D.eq
  return AffinePlaneCurve(F)
end

####################################################
end
using .AffinePlaneCurveModule
