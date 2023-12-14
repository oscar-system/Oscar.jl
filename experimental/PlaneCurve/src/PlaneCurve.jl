module PlaneCurveModule
using Oscar
import Base.==

export AffinePlaneCurve
export Point
export ProjPlaneCurve
export ProjectivePlaneCurve
export curve_components
export defining_equation
export degree
export hash
export ideal_point
export is_irreducible
export is_reduced
export jacobi_ideal
export reduction
export ring
export union

################################################################################

abstract type PlaneCurve{S <: RingElem} end
abstract type ProjectivePlaneCurve{S} <: PlaneCurve{S} end

################################################################################
# Point (not specific to curves).
################################################################################
# Point (first attempt).
# The point is described by its coordinates.
# ambient_dim gives the dimension of the space to which belongs the point.
@doc raw"""
    Point(coordinates::Vector{S}) where {S <: FieldElem}

Return the point with the given coordinates.

# Examples
```julia
julia> P = Point([QQ(1), QQ(2), QQ(2)])
Point with coordinates QQFieldElem[1, 2, 2]
```
"""
mutable struct Point{S <: FieldElem}
  coord::Vector{S}
  ambient_dim::Int
  function Point(coordinates::Vector{S}) where {S <: FieldElem}
    r = new{S}()
    r.coord = coordinates
    r.ambient_dim = length(coordinates)
    return r
  end
end

function Base.show(io::IO, P::Point)
  print(io, "Point with coordinates ", P.coord)
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

@doc raw"""
    ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem

Return the maximal ideal associated to the point `P` in the ring `R`.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> P = Point([QQ(2), QQ(1)])
Point with coordinates QQFieldElem[2, 1]

julia> Oscar.ideal_point(R, P)
ideal(x - 2, y - 1)
```
"""
function ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem
  V = gens(R)
  return ideal(R, [V[i] - P.coord[i] for i in 1:length(V)])
end

################################################################################
# Structure of Affine Plane Curves and Projective Plane Curves
################################################################################
@doc raw"""
    AffinePlaneCurve{S}(eq::Oscar.MPolyRingElem{S}) where S <: FieldElem

Return the Affine Plane Curve defined by the polynomial in two variables `eq`.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> F = y^3*x^6 - y^6*x^2
x^6*y^3 - x^2*y^6

julia> C = AffinePlaneCurve(F)
Affine plane curve defined by x^6*y^3 - x^2*y^6
```
"""
mutable struct AffinePlaneCurve{S} <: PlaneCurve{S}
  eq::Oscar.MPolyRingElem{S}                # Equation of the curve (polynomial in two variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{AffinePlaneCurve{S}, Int}
  function AffinePlaneCurve{S}(eq::Oscar.MPolyRingElem{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 2 || error("The defining equation must belong to a ring with two variables")
    !is_constant(eq) || error("The defining equation must be non constant")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{AffinePlaneCurve{S}, Int}())
  end
end

AffinePlaneCurve(eq::Oscar.MPolyRingElem{S}) where {S <: FieldElem} = AffinePlaneCurve{S}(eq)

function Base.show(io::IO, C::AffinePlaneCurve)
  if !get(io, :compact, false)
     print(io, "Affine plane curve defined by ", C.eq)
  else
     print(io, C.eq)
  end
end

################################################################################
@doc raw"""
    ProjPlaneCurve{S}(eq::Oscar.MPolyDecRingElem{S}) where {S <: FieldElem}

Return the Projective Plane Curve defined by the homogeneous polynomial in three variables `eq`.

# Examples
```julia
julia> R, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> T, _ = grade(R)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> F = T(y^3*x^6 - y^6*x^2*z)
x^6*y^3 - x^2*y^6*z

julia> ProjPlaneCurve(F)
Projective plane curve defined by x^6*y^3 - x^2*y^6*z
```
"""
mutable struct ProjPlaneCurve{S} <: ProjectivePlaneCurve{S}
  eq::Oscar.MPolyDecRingElem{S}            # Equation of the curve (polynomial in three variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{ProjPlaneCurve{S}, Int}
  function ProjPlaneCurve{S}(eq::Oscar.MPolyDecRingElem{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
    !is_constant(eq) || error("The defining equation must be non constant")
    is_homogeneous(eq) || error("The defining equation is not homogeneous")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{ProjPlaneCurve{S}, Int}())
  end
end

ProjPlaneCurve(eq::Oscar.MPolyDecRingElem{S}) where {S <: FieldElem} = ProjPlaneCurve{S}(eq)

@doc raw"""
    ProjPlaneCurve(f::MPolyRingElem{T}) where {T <: FieldElem}

Given a homogeneous polynomial `f` in three variables with coefficients in a field,
create the projective plane curve defined by `f`.

# Examples
```julia
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> C = ProjPlaneCurve(z*x^2-y^3)
Projective plane curve defined by x^2*z - y^3
```
"""
function ProjPlaneCurve(eq::Oscar.MPolyRingElem{S}) where {S <: FieldElem}
  R, _ = grade(parent(eq))
  return ProjPlaneCurve{S}(R(eq))
end

function Base.show(io::IO, C::ProjPlaneCurve)
  if !get(io, :compact, false)
     print(io, "Projective plane curve defined by ", C.eq)
  else
     print(io, C.eq)
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
  F = 1//AbstractAlgebra.leading_coefficient(C.eq)*C.eq
  return hash(F, h)
end

################################################################################
# check for equality of curves (the equations are equal up to multiplication by
# a non zero constant).

function ==(C::PlaneCurve{S}, D::PlaneCurve{S}) where S <: FieldElem
  F = defining_equation(C)
  G = defining_equation(D)
  return degree(C) == degree(D) && F*(AbstractAlgebra.leading_coefficient(G)//AbstractAlgebra.leading_coefficient(F)) == G
end

################################################################################
@doc raw"""
    in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.
"""
function Base.in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem
  return iszero(evaluate(C.eq, P.coord))
end

@doc raw"""
    in(P::Oscar.Geometry.ProjSpcElem{S}, C::ProjectivePlaneCurve{S}) where S <: FieldElem

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.
"""
function Base.in(P::Oscar.Geometry.ProjSpcElem{S}, C::ProjectivePlaneCurve{S}) where S <: RingElem
  return iszero(evaluate(C.eq, P.v))
end

################################################################################
# Compute the degree of the equation of the curve if not already known,
# and show it.

@doc raw"""
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

@doc raw"""
    jacobi_ideal(C::PlaneCurve)

Return the Jacobian ideal of the defining polynomial of `C`.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.jacobi_ideal(C)
ideal(6*x^5*y^3 - 2*x*y^6, 3*x^6*y^2 - 6*x^2*y^5)
```
"""
function Oscar.jacobi_ideal(C::PlaneCurve)
 return jacobi_ideal(C.eq)
end

################################################################################
# Components of the curve

@doc raw"""
    curve_components(C::PlaneCurve{S}) where S <: FieldElem

Return a dictionary containing the irreducible components of `C` and their multiplicity.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.curve_components(C)
Dict{AffinePlaneCurve{QQFieldElem}, Int64} with 3 entries:
  y         => 3
  x         => 2
  x^4 - y^3 => 1
```
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

@doc raw"""
    is_irreducible(C::PlaneCurve{S}) where S <: FieldElem

Return `true` if `C` is irreducible, and `false` otherwise.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> Oscar.is_irreducible(C)
true

julia> D = AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.is_irreducible(D)
false
```
"""
function Oscar.is_irreducible(C::PlaneCurve{S}) where S <: FieldElem
   return is_irreducible(defining_equation(C))
end

################################################################################
# Check reducedness by computing a factorization

@doc raw"""
    is_reduced(C::PlaneCurve{S}) where S <: FieldElem

Return `true` if `C` is reduced, and `false` otherwise.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> Oscar.is_reduced(C)
true

julia> D = AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.is_reduced(D)
false
```
"""
function Oscar.is_reduced(C::PlaneCurve{S}) where S <: FieldElem
  if isempty(C.components)
     L = factor_squarefree(defining_equation(C))
     return all(isone, values(L.fac))
  else
     return all(isone, values(C.components))
  end
end

################################################################################
# Compute squarefree equation.

@doc raw"""
    reduction(C::AffinePlaneCurve{S}) where S <: FieldElem
    reduction(C::ProjPlaneCurve{S}) where S <: FieldElem

Return the plane curve defined by the squarefree part of the equation of `C`.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.reduction(C)
Affine plane curve defined by x^5*y - x*y^4
```
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

@doc raw"""
    union(C::T, D::T) where T <: PlaneCurve

Return the union of `C` and `D` (with multiplicity).

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> D = AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> union(C, D)
Affine plane curve defined by -x^9*y^3 + x^7*y^3 + x^6*y^5 + x^5*y^6 - x^3*y^6 - x^2*y^8
```
"""
Base.union(C::T, D::T) where T <: PlaneCurve = T(C.eq*D.eq)

################################################################################
# Ring associated to a curve

@doc raw"""
    ring(C::PlaneCurve)

Return the coordinate ring of the curve `C`.

# Examples
```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> Oscar.ring(C)
(Quotient of multivariate polynomial ring by ideal(-x^3 + x + y^2), Map: multivariate polynomial ring -> quotient of multivariate polynomial ring)
```
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
include("EllCurveZnZ.jl")
include("ProjCurve.jl")
include("ParaPlaneCurves.jl")

################################################################################
end
#=
using .PlaneCurveModule

export Point, ideal_point, AffinePlaneCurve, ProjPlaneCurve, hash, degree,
       jacobi_ideal, curve_components, is_irreducible, is_reduced, reduction,
       union, defining_equation, ring, ProjectivePlaneCurve

export is_smooth, tangent, common_components, curve_intersect,
       curve_singular_locus, is_smooth_curve, multiplicity,
       tangent_lines, intersection_multiplicity, aretransverse,
       arithmetic_genus, geometric_genus

export ProjCurve, defining_ideal, curve_components, reduction, is_irreducible,
       jacobi_ideal
       
export parametrization_plane_curve, adjoint_ideal, rational_point_conic,
       parametrization_conic, map_to_rational_normal_curve,
       rat_normal_curve_anticanonical_map, rat_normal_curve_It_Proj_Odd,
       rat_normal_curve_It_Proj_Even, invert_birational_map

=#
