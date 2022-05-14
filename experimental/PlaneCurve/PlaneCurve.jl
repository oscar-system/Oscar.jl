module PlaneCurveModule
using Oscar, Markdown
import Base.==

export Point, ideal_point, AffinePlaneCurve, ProjPlaneCurve, hash, degree,
       jacobi_ideal, curve_components, is_irreducible, isreduced, reduction,
       union, defining_equation, ring, ProjectivePlaneCurve

################################################################################

abstract type PlaneCurve{S <: RingElem} end
abstract type ProjectivePlaneCurve{S} <: PlaneCurve{S} end

################################################################################
# Point (not specific to curves).
################################################################################
# Point (first attempt).
# The point is described by its coordinates.
# ambient_dim gives the dimension of the space to which belongs the point.
@doc Markdown.doc"""
    Point(coordinates::Vector{S}) where {S <: FieldElem}

Return the point with the given coordinates.

# Examples
```jldoctest
julia> P = Oscar.Point([QQ(1), QQ(2), QQ(2)])
Point with coordinates fmpq[1, 2, 2]
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

@doc Markdown.doc"""
    ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem

Return the maximal ideal associated to the point `P` in the ring `R`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> P = Oscar.Point([QQ(2), QQ(1)])
Point with coordinates fmpq[2, 1]

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
@doc Markdown.doc"""
    AffinePlaneCurve{S}(eq::Oscar.MPolyElem{S}) where S <: FieldElem

Return the Affine Plane Curve defined by the polynomial in two variables `eq`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = y^3*x^6 - y^6*x^2
x^6*y^3 - x^2*y^6

julia> C = Oscar.AffinePlaneCurve(F)
Affine plane curve defined by x^6*y^3 - x^2*y^6
```
"""
mutable struct AffinePlaneCurve{S} <: PlaneCurve{S}
  eq::Oscar.MPolyElem{S}                # Equation of the curve (polynomial in two variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{AffinePlaneCurve{S}, Int}
  function AffinePlaneCurve{S}(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 2 || error("The defining equation must belong to a ring with two variables")
    !is_constant(eq) || error("The defining equation must be non constant")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{AffinePlaneCurve{S}, Int}())
  end
end

AffinePlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem} = AffinePlaneCurve{S}(eq)

function Base.show(io::IO, C::AffinePlaneCurve)
  if !get(io, :compact, false)
     print(io, "Affine plane curve defined by ", C.eq)
  else
     print(io, C.eq)
  end
end

################################################################################
@doc Markdown.doc"""
    ProjPlaneCurve{S}(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}

Return the Projective Plane Curve defined by the homogeneous polynomial in three variables `eq`.

# Examples
```jldoctest
julia> R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(R)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> F = T(y^3*x^6 - y^6*x^2*z)
x^6*y^3 - x^2*y^6*z

julia> Oscar.ProjPlaneCurve(F)
Projective plane curve defined by x^6*y^3 - x^2*y^6*z
```
"""
mutable struct ProjPlaneCurve{S} <: ProjectivePlaneCurve{S}
  eq::Oscar.MPolyElem_dec{S}            # Equation of the curve (polynomial in three variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{ProjPlaneCurve{S}, Int}
  function ProjPlaneCurve{S}(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
    !is_constant(eq) || error("The defining equation must be non constant")
    is_homogeneous(eq) || error("The defining equation is not homogeneous")
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{ProjPlaneCurve{S}, Int}())
  end
end

ProjPlaneCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem} = ProjPlaneCurve{S}(eq)

@doc Markdown.doc"""
    ProjPlaneCurve(f::MPolyElem{T}) where {T <: FieldElem}

Given a homogeneous polynomial `f` in three variables with coefficients in a field,
create the projective plane curve defined by `f`.

# Examples
```jldoctest
julia> R, (x,y,z) = GradedPolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = ProjPlaneCurve(z*x^2-y^3)
Projective plane curve defined by x^2*z - y^3
```
"""
function ProjPlaneCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
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
  F = 1//leading_coefficient(C.eq)*C.eq
  return hash(F, h)
end

################################################################################
# check for equality of curves (the equations are equal up to multiplication by
# a non zero constant).

function ==(C::PlaneCurve{S}, D::PlaneCurve{S}) where S <: FieldElem
  F = defining_equation(C)
  G = defining_equation(D)
  return degree(C) == degree(D) && F*(leading_coefficient(G)//leading_coefficient(F)) == G
end

################################################################################
@doc Markdown.doc"""
    in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.
"""
function Base.in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem
  return iszero(evaluate(C.eq, P.coord))
end

@doc Markdown.doc"""
    in(P::Oscar.Geometry.ProjSpcElem{S}, C::ProjectivePlaneCurve{S}) where S <: FieldElem

Return `true` if the point `P` is on the curve `C`, and `false` otherwise.
"""
function Base.in(P::Oscar.Geometry.ProjSpcElem{S}, C::ProjectivePlaneCurve{S}) where S <: RingElem
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

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^3*x^6 - y^6*x^2)
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

@doc Markdown.doc"""
    curve_components(C::PlaneCurve{S}) where S <: FieldElem

Return a dictionary containing the irreducible components of `C` and their multiplicity.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.curve_components(C)
Dict{AffinePlaneCurve{fmpq}, Int64} with 3 entries:
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

@doc Markdown.doc"""
    is_irreducible(C::PlaneCurve{S}) where S <: FieldElem

Return `true` if `C` is irreducible, and `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> Oscar.is_irreducible(C)
true

julia> D = Oscar.AffinePlaneCurve(y^3*x^6 - y^6*x^2)
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

@doc Markdown.doc"""
    isreduced(C::PlaneCurve{S}) where S <: FieldElem

Return `true` if `C` is reduced, and `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> Oscar.isreduced(C)
true

julia> D = Oscar.AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> Oscar.isreduced(D)
false
```
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

@doc Markdown.doc"""
    reduction(C::AffinePlaneCurve{S}) where S <: FieldElem
    reduction(C::ProjPlaneCurve{S}) where S <: FieldElem

Return the plane curve defined by the squarefree part of the equation of `C`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^3*x^6 - y^6*x^2)
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

@doc Markdown.doc"""
    union(C::T, D::T) where T <: PlaneCurve

Return the union of `C` and `D` (with multiplicity).

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> D = Oscar.AffinePlaneCurve(y^3*x^6 - y^6*x^2)
Affine plane curve defined by x^6*y^3 - x^2*y^6

julia> union(C, D)
Affine plane curve defined by -x^9*y^3 + x^7*y^3 + x^6*y^5 + x^5*y^6 - x^3*y^6 - x^2*y^8
```
"""
Base.union(C::T, D::T) where T <: PlaneCurve = T(C.eq*D.eq)

################################################################################
# Ring associated to a curve

@doc Markdown.doc"""
    ring(C::PlaneCurve)

Return the coordinate ring of the curve `C`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> C = Oscar.AffinePlaneCurve(y^2+x-x^3)
Affine plane curve defined by -x^3 + x + y^2

julia> Oscar.ring(C)
(Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(-x^3 + x + y^2), Map from
Multivariate Polynomial Ring in x, y over Rational Field to Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(-x^3 + x + y^2) defined by a julia-function with inverse)
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
using .PlaneCurveModule

export Point, ideal_point, AffinePlaneCurve, ProjPlaneCurve, hash, degree,
       jacobi_ideal, curve_components, is_irreducible, isreduced, reduction,
       union, defining_equation, ring, ProjectivePlaneCurve

export issmooth, tangent, common_components, curve_intersect,
       curve_singular_locus, issmooth_curve, multiplicity,
       tangent_lines, intersection_multiplicity, aretransverse,
       arithmetic_genus, geometric_genus

export ProjCurve, defining_ideal, curve_components, reduction, is_irreducible,
       jacobi_ideal
       
export parametrization_plane_curve, adjoint_ideal, rational_point_conic,
       parametrization_conic, map_to_rational_normal_curve,
       rat_normal_curve_anticanonical_map, rat_normal_curve_It_Proj_Odd,
       rat_normal_curve_It_Proj_Even, invert_birational_map

