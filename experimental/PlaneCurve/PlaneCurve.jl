module PlaneCurveModule
using Oscar, Markdown
import Base.==

export Point, ideal_point, AffinePlaneCurve, ProjPlaneCurve, hash, degree,
       jacobi_ideal, curve_components, isirreducible, isreduced, reduction,
       union, defining_equation, ring

################################################################################

abstract type PlaneCurve{S <: Union{RingElem, FieldElem}} end
abstract type ProjectivePlaneCurve{S} <: PlaneCurve{S} end

################################################################################
# Point (not specific to curves).
################################################################################
# Point (first attempt).
# The point is described by its coordinates.
# ambient_dim gives the dimension of the space to which belongs the point.

mutable struct Point{S <: FieldElem}
  coord::Vector{S}
  ambient_dim::Int
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
# check for equality of points (the coordinates are equal).

function ==(P::Point, Q::Point)
  return P.coord == Q.coord
end

################################################################################

function Base.hash(P::Point, h::UInt)
  return hash(P.coord, h)
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

mutable struct AffinePlaneCurve{S} <: PlaneCurve{S}
  eq::Oscar.MPolyElem{S}                # Equation of the curve (polynomial in two variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{AffinePlaneCurve{S}, Int}
  function AffinePlaneCurve{S}(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 2 || error("The defining equation must belong to a ring with two variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{AffinePlaneCurve{S}, Int}())
  end
end

AffinePlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem} = AffinePlaneCurve{S}(eq)

function Base.show(io::IO, C::AffinePlaneCurve)
  if !get(io, :compact, false)
     println(io, "Affine plane curve defined by ", C.eq)
  else
     println(io, C.eq)
  end
end

################################################################################

mutable struct ProjPlaneCurve{S} <: ProjectivePlaneCurve{S}
  eq::Oscar.MPolyElem_dec{S}            # Equation of the curve (polynomial in three variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{ProjPlaneCurve{S}, Int}
  function ProjPlaneCurve{S}(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    ishomogenous(eq) || error("The defining equation is not homogeneous")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{ProjPlaneCurve{S}, Int}())
  end
end

ProjPlaneCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem} = ProjPlaneCurve{S}(eq)
function ProjPlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
  R = grade(parent(eq))
  return ProjPlaneCurve{S}(R(eq))
end

function Base.show(io::IO, C::ProjPlaneCurve)
  if !get(io, :compact, false)
     println(io, "Projective plane curve defined by ", C.eq)
  else
     println(io, C.eq)
  end
end

################################################################################
################################################################################
# Plane Curves related functions.
################################################################################

defining_equation(C::PlaneCurve) = C.eq
defining_equation(C::ProjectivePlaneCurve) = C.eq.f


Oscar.dim(::PlaneCurve) = 1 # since C is a plane curve, the dimension is always 1


################################################################################
# hash function

function Base.hash(C::PlaneCurve, h::UInt)
  F = 1//leading_coefficient(C.eq)*C.eq
  return hash(F, h)
end

################################################################################
# check for equality of curves (the equations are equal up to multiplication by
# a non zero constant).

function ==(C::PlaneCurve, D::PlaneCurve)
  F = defining_equation(C)
  G = defining_equation(D)
  return degree(C) == degree(D) && F*(leading_coefficient(G)//leading_coefficient(F)) == G
end

################################################################################
@doc Markdown.doc"""
    in(P::Point{S}, C::AffinePlaneCurve{S})

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.
"""
function Base.in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem
  return iszero(evaluate(C.eq, P.coord))
end

@doc Markdown.doc"""
    in(P::Oscar.Geometry.ProjSpcElem{S}, C::ProjectivePlaneCurve{S}) where S <: FieldElem

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.
"""
function Base.in(P::Oscar.Geometry.ProjSpcElem{S}, C::ProjectivePlaneCurve{S}) where S <: FieldElem
  return iszero(evaluate(C.eq, P.v))
end

################################################################################
# Compute the degree of the equation of the curve if not already known,
# and show it.

@doc Markdown.doc"""
    degree(C::PlaneCurve)

Return the degree of the defining polynomial of `C`.
"""
function Oscar.degree(C::PlaneCurve)
  if C.degree < 0
    C.degree = total_degree(defining_equation(C))
  end
  return C.degree
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
# Components of the curve

@doc Markdown.doc"""
    curve_components(C::PlaneCurve{S}) where S <: FieldElem

Return a dictionary containing the irreducible components of `C` and their multiplicity.
"""
function curve_components(C::PlaneCurve{S}) where S <: FieldElem
  if isempty(C.components)
    T = typeof(C)
    D = factor(C.eq)
    C.components = Dict(T(x) => D.fac[x] for x in keys(D.fac))
  end
  return C.components
end

################################################################################
# Check irreducibility.
# TODO: change for a direct irreducibility check when available.

@doc Markdown.doc"""
    isirreducible(C::PlaneCurve{S}) where S <: FieldElem

Return `true` if `C` is irreducible, and `false` otherwise.
"""
function Oscar.isirreducible(C::PlaneCurve{S}) where S <: FieldElem
   comp = curve_components(C)
   return length(comp) == 1 && all(isone, values(comp))
end

################################################################################
# Check reducedness by computing a factorization

@doc Markdown.doc"""
    isreduced(C::PlaneCurve{S}) where S <: FieldElem

Return `true` if `C` is reduced, and `false` otherwise.
"""
function Oscar.isreduced(C::PlaneCurve{S}) where S <: FieldElem
  if isempty(C.components)
     L = factor_squarefree(defining_equation(C))
     return all(isone, values(L.fac))
  else
     return all(isone, values(C.components))
  end
end

################################################################################
# Compute squarefree equation.
# TODO: change for a direct squarefree computation when available.

@doc Markdown.doc"""
    reduction(C::PlaneCurve{S}) where S <: FieldElem

Return the plane curve defined by the squarefree part of the equation of `C`.
"""
function reduction(C::AffinePlaneCurve{S}) where S <: FieldElem
  if isempty(C.components)
     L = factor_squarefree(C.eq)
     F = prod(f -> f, keys(L.fac))
     rC = AffinePlaneCurve(F)
     return rC
  else
     F = prod(D -> D.eq, keys(C.components))
     rC = AffinePlaneCurve(F)
     rC.components = Dict(AffinePlaneCurve(D.eq) => 1 for D in keys(C.components))
     return rC
  end
end

function reduction(C::ProjPlaneCurve{S}) where S <: FieldElem
  comp = curve_components(C)
  F = prod(D -> D.eq, keys(comp))
  rC = ProjPlaneCurve(F)
  rC.components = Dict(ProjPlaneCurve(D.eq) => 1 for D in keys(comp))
  return rC
end

################################################################################
# Union of two plane curves of the same type (with multiplicity)

@doc Markdown.doc"""
   union(C::T, D::T) where T <: PlaneCurve

Return the union of `C` and `D` (with multiplicity).
"""
Base.union(C::T, D::T) where T <: PlaneCurve = T(C.eq*D.eq)

################################################################################
# Ring associated to a curve

@doc Markdown.doc"""
    ring(C::PlaneCurve)

Return the coordinate ring of the curve `C`.
"""
function ring(C::PlaneCurve)
  F = C.eq
  S = parent(F)
  return quo(S, ideal(S, [F]))
end

################################################################################

include("AffinePlaneCurve.jl")
include("ProjPlaneCurve.jl")
include("DivisorCurve.jl")
include("ProjEllipticCurve.jl")

################################################################################
end
using .PlaneCurveModule
