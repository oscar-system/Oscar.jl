module PlaneCurveModule
using Oscar, Markdown
import Base.:(==)

export Point, ideal_point, AffinePlaneCurve, ProjPlaneCurve, hash, degree,
       jacobi_ideal, curve_components, isirreducible, isreduced, reduction,
       union

################################################################################

abstract type PlaneCurve end

################################################################################
# Point (not specific to curves).
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
# Associate a maximal ideal to a point in a given ring (not specific to curves)

@doc Markdown.doc"""
    ideal_point(P::Point{S}, R::MPolyRing{S}) where S <: FieldElem

Return the maximal ideal associated to the point `P` in the ring `R`.
"""
function ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem
  V = gens(R)
  return ideal(R, [V[i] - P.coord[i] for i in 1:length(V)])
end

################################################################################
# Structure of Affine Plane Curves and Projective Plane Curves
################################################################################

mutable struct AffinePlaneCurve{S <: FieldElem} <: PlaneCurve
  eq::Oscar.MPolyElem{S}                  # Equation of the curve (polynomial in two variables)
  degree::Int64                           # degree of the equation of the curve
  dimension::Int64                        # dimension of the curve (as a variety)
  components::Dict
  function AffinePlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    if nvars(parent(eq)) != 2
       error("The defining equation must belong to a ring with two variables")
    elseif isconstant(eq)
       error("The defining equation must be non constant")
    end
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            1,                   # since C is a plane curve, the dimension is always 1
            Dict{AffinePlaneCurve{S}, Int64}())
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

mutable struct ProjPlaneCurve{S <: FieldElem} <: PlaneCurve
  eq::Oscar.MPolyElem_dec{S}              # Equation of the curve (polynomial in three variables)
  degree::Int64                           # degree of the equation of the curve
  dimension::Int64                        # dimension of the curve (as a variety)
  components::Dict{ProjPlaneCurve{S}, Int64}
  function ProjPlaneCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    if nvars(parent(eq)) != 3
       error("The defining equation must belong to a ring with three variables")
    elseif isconstant(eq)
       error("The defining equation must be non constant")
    elseif !ishomogenous(eq)
       error("The defining equation is not homogeneous")
    end
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            1,                   # since C is a plane curve, the dimension is always 1
            Dict{ProjPlaneCurve{S}, Int64}())
  end
  function ProjPlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
     R = grade(parent(eq))
     return ProjPlaneCurve(R(eq))
  end
end

function Base.show(io::IO, C::ProjPlaneCurve)
  if !get(io, :compact, false)
     println(io, "Projective plane curve defined by ", C.eq)
  else
     println(io, C.eq)
  end
end

################################################################################
# Plane Curves related functions.
################################################################################
# hash function

function Base.hash(C::PlaneCurve, h::UInt)
  F = 1//lc(C.eq)*C.eq
  return hash(F, h)
end

################################################################################
# check for equality of curves (the equations are equal up to multiplication by
# a non zero constant).

function ==(C::PlaneCurve, D::PlaneCurve)
  F = C.eq
  G = D.eq
  return degree(C) == degree(D) && F*(lc(G)//lc(F)) == G
end

################################################################################
# Compute the degree of the equation of the curve if not already known,
# and show it.

@doc Markdown.doc"""
    degree(C::PlaneCurve)

Return the degree of the defining polynomial of `C`.
"""
function Oscar.degree(C::PlaneCurve)
  if C.degree < 0 && typeof(C) == AffinePlaneCurve{typeof(lc(C.eq))}
    C.degree = total_degree(C.eq)
  elseif C.degree < 0 && typeof(C) == ProjPlaneCurve{typeof(lc(C.eq))}
    C.degree = total_degree(C.eq.f)
  end
  return(C.degree)
end

################################################################################

@doc Markdown.doc"""
    jacobi_ideal(C::PlaneCurve)

Return the Jacobian ideal of the defining polynomial of `C`.
"""
function Oscar.jacobi_ideal(C::PlaneCurve)
 return jacobi_ideal(C.eq)
end

################################################################################
# helping function: computes the irreducible components of the curve

function compo(C::PlaneCurve)
  if typeof(C) == AffinePlaneCurve{typeof(lc(C.eq))}
    D = factor(C.eq)
    C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  elseif typeof(C) == ProjPlaneCurve{typeof(lc(C.eq))}
    D = factor(C.eq.f)
    C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
    curve_components(C::PlaneCurve)

Return a dictionary containing the irreducible components of `C` and their multiplicity.
"""
function curve_components(C::PlaneCurve)
  if isempty(C.components)
    compo(C)
  end
  return C.components
end

################################################################################
# Check irreducibility.
# TODO: change for a direct irreducibility check when available.

@doc Markdown.doc"""
    isirreducible(C::PlaneCurve)

Return `true` if `C` is irreducible, and `false` otherwise.
"""
function Oscar.isirreducible(C::PlaneCurve)
  if isempty(C.components)
     compo(C)
  end
  if length(C.components) != 1
     return false
  else
     return all(isone, values(C.components))
  end
end

################################################################################
# Check reducedness by computing a factorization
# TODO: change for a direct squarefree check when available.

@doc Markdown.doc"""
    isreduced(C::PlaneCurve)

Return `true` if `C` is reduced, and `false` otherwise.
"""
function Oscar.isreduced(C::PlaneCurve)
  if isempty(C.components)
    compo(C)
  end
  return all(isone, values(C.components))
end

################################################################################
# Compute squarefree equation.
# TODO: change for a direct squarefree computation when available.

@doc Markdown.doc"""
    reduction(C::PlaneCurve)

Return the plane curve defined by the squarefree part of the equation of `C`.
"""
function reduction(C::PlaneCurve)
  if isempty(C.components)
     compo(C)
  end
  F = prod(D -> D.eq, keys(C.components))
  if typeof(C) == AffinePlaneCurve{typeof(lc(C.eq))}
    rC = AffinePlaneCurve(F)
    rC.components = Dict(AffinePlaneCurve(D.eq) => 1 for D in keys(C.components))
  elseif typeof(C) == ProjPlaneCurve{typeof(lc(C.eq))}
    rC = ProjPlaneCurve(F)
    rC.components = Dict(ProjPlaneCurve(D.eq) => 1 for D in keys(C.components))
  end
  return rC
end

################################################################################
# Union of two plane curves of the same type (with multiplicity)

@doc Markdown.doc"""
   union(C::PlaneCurve{S}, D::PlaneCurve{S}) where S <: FieldElem

Return the union of `C` and `D` (with multiplicity).
"""
function Base.union(C::PlaneCurve, D::PlaneCurve)
  typeof(C) == typeof(D) || error("The curves do not have the same type.")
  F = C.eq*D.eq
  if typeof(C) == ProjPlaneCurve{typeof(lc(C.eq))}
     return ProjPlaneCurve(F)
  elseif typeof(C) == AffinePlaneCurve{typeof(lc(C.eq))}
     return AffinePlaneCurve(F)
  end
end

################################################################################

include("AffinePlaneCurve.jl")
include("ProjPlaneCurve.jl")

################################################################################
end
using .PlaneCurveModule
