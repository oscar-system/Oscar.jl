module PlaneCurveModule
using Oscar, Markdown
import Base.:(==)

export factor, gcd, div, Point, ideal_point, AffinePlaneCurve, ProjPlaneCurve,
       hash, degree, jacobi_ideal, curve_components, isirreducible, isreduced,
       reduction, union, check_on_curve

################################################################################

abstract type PlaneCurve end

################################################################################
# TODO: Move these three functions to the file mpoly-graded.jl
################################################################################

function Oscar.factor(x::Oscar.MPolyElem_dec)
  R = parent(x)
  D = Dict{elem_type(R), Int64}()
  F = factor(x.f)
  n=length(F.fac)
  #if n == 1
  #  return Fac(R(F.unit), D)
  #else
    for i in keys(F.fac)
     push!(D, R(i) => Int64(F[i]))
    end
  return Fac(R(F.unit), D)
  #end
end


function Oscar.gcd(x::Oscar.MPolyElem_dec, y::Oscar.MPolyElem_dec)
  R = parent(x)
  return R(gcd(x.f, y.f))
end

function Oscar.div(x::Oscar.MPolyElem_dec, y::Oscar.MPolyElem_dec)
  R = parent(x)
  return R(div(x.f, y.f))
end

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

mutable struct AffinePlaneCurve{S <: FieldElem} <: PlaneCurve
  eq::Oscar.MPolyElem{S}                # Equation of the curve (polynomial in two variables)
  degree::Int                           # degree of the equation of the curve
  dimension::Int                        # dimension of the curve (as a variety)
  components::Dict{AffinePlaneCurve{S}, Int}
  function AffinePlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 2 || error("The defining equation must belong to a ring with two variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            1,                   # since C is a plane curve, the dimension is always 1
            Dict{AffinePlaneCurve{S}, Int}())
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
  eq::Oscar.MPolyElem_dec{S}            # Equation of the curve (polynomial in three variables)
  degree::Int                           # degree of the equation of the curve
  dimension::Int                        # dimension of the curve (as a variety)
  components::Dict{ProjPlaneCurve{S}, Int}
  function ProjPlaneCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    ishomogenous(eq) || error("The defining equation is not homogeneous")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            1,                   # since C is a plane curve, the dimension is always 1
            Dict{ProjPlaneCurve{S}, Int}())
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
################################################################################
# Plane Curves related functions.
################################################################################

defining_equation(C::AffinePlaneCurve) = C.eq
defining_equation(C::ProjPlaneCurve) = C.eq.f



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
  F = defining_equation(C)
  G = defining_equation(D)
  return degree(C) == degree(D) && F*(lc(G)//lc(F)) == G
end

################################################################################

function check_on_curve(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  return iszero(evaluate(C.eq, P.coord))
end

function check_on_curve(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
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
# helping function: computes the irreducible components of the curve

function compo(C::AffinePlaneCurve)
  D = factor(C.eq)
  C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
end

function compo(C::ProjPlaneCurve)
  D = factor(C.eq)
  C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
end

################################################################################
# Helping function

function _assure_has_components(C::PlaneCurve)
  if isempty(C.components)
     compo(C)
  end
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
    curve_components(C::PlaneCurve)

Return a dictionary containing the irreducible components of `C` and their multiplicity.
"""
function curve_components(C::PlaneCurve)
  _assure_has_components(C)
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
  _assure_has_components(C)
   return length(C.components) == 1 && all(isone, values(C.components))
end

################################################################################
# Check reducedness by computing a factorization

@doc Markdown.doc"""
    isreduced(C::PlaneCurve)

Return `true` if `C` is reduced, and `false` otherwise.
"""
function Oscar.isreduced(C::PlaneCurve)
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
    reduction(C::PlaneCurve)

Return the plane curve defined by the squarefree part of the equation of `C`.
"""
function reduction(C::AffinePlaneCurve)
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

function reduction(C::ProjPlaneCurve)
  _assure_has_components(C)
  F = prod(D -> D.eq, keys(C.components))
  rC = ProjPlaneCurve(F)
  rC.components = Dict(ProjPlaneCurve(D.eq) => 1 for D in keys(C.components))
  return rC
end

################################################################################
# Union of two plane curves of the same type (with multiplicity)

@doc Markdown.doc"""
   union(C::PlaneCurve{S}, D::PlaneCurve{S}) where S <: FieldElem

Return the union of `C` and `D` (with multiplicity).
"""
function Base.union(C::T, D::T)  where T <: AffinePlaneCurve
     return AffinePlaneCurve(C.eq*D.eq)
end
function Base.union(C::T, D::T)  where T <: ProjPlaneCurve
     return ProjPlaneCurve(C.eq*D.eq)
end

################################################################################

include("AffinePlaneCurve.jl")
include("ProjPlaneCurve.jl")
include("DivisorCurve.jl")

################################################################################
end
using .PlaneCurveModule
