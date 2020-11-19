module ProjPlaneCurveModule

using Oscar, Markdown

include("AffinePlaneCurve.jl")
using .AffinePlaneCurveModule

import Base.:(==)

export ProjPlaneCurve, degree, jacobi_ideal, issmooth, tangent,
       curve_components, isirreducible, isreduced, common_component, reduction,
       dehomogenization, curve_intersect, curve_singular_locus, issmooth_curve

################################################################################
#
# We follow the lecture notes by Andreas Gathmann:
# https://www.mathematik.uni-kl.de/~gathmann/en/curves.php
#
################################################################################
#

mutable struct ProjPlaneCurve{S <: FieldElem}
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
# hash function

function Base.hash(C::ProjPlaneCurve, h::UInt)
  F = 1//lc(C.eq)*C.eq
  return hash(F, h)
end

################################################################################

################################################################################
# Compute the degree of the equation of the curve if not already known,
# and show it.

@doc Markdown.doc"""
    degree(C::ProjPlaneCurve)

Return the degree of the defining polynomial of `C`.
"""
function Oscar.degree(C::ProjPlaneCurve)
   if C.degree < 0
      C.degree = total_degree(C.eq.f)
   end
   return C.degree
end

################################################################################
# Compute the Jacobian Ideal of C

@doc Markdown.doc"""
    jacobi_ideal(C::ProjPlaneCurve)

Return the Jacobian ideal of the defining polynomial of `C`.
"""
function Oscar.jacobi_ideal(C::ProjPlaneCurve)
 return jacobi_ideal(C.eq)
end

################################################################################
# To check if a point is smooth: return true if the point is a smooth point on
# the curve C, and false  if it is singular, and an error if it is not on the
# curve.

@doc Markdown.doc"""
    issmooth(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Throw an error if `P` is not a point of `C`, return `false` if `P` is a singular point of `C`, and `true` if `P` is a smooth point of `C`.
"""
function Oscar.issmooth(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  dim(P.parent) == 2 || error("The point needs to be in a projective two dimensional space")
  iszero(evaluate(C.eq, P.v)) || error("The point is not on the curve defined by ", C.eq)
  J = jacobi_ideal(C)
  L = gens(J)
  a = evaluate(L[1], P.v)
  b = evaluate(L[2], P.v)
  c = evaluate(L[3], P.v)
  return !(iszero(a) && iszero(b) && iszero(c))
end

################################################################################
# If P is a smooth point of a projective plane curve C, compute the tangent of C
# at the point P.

@doc Markdown.doc"""
    tangent(C::ProjPlaneCurve{S}, P::Point{S}) where S <: FieldElem

Return the tangent of `C` at `P` when `P` is a smooth point of `C`, and throw an error otherwise.
"""
function tangent(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  dim(P.parent) == 2 || error("The point needs to be in a projective two dimensional space")
  iszero(evaluate(C.eq, P.v)) || error("The point is not on the curve")
  J = jacobi_ideal(C)
  L = gens(J)
  a = evaluate(L[1], P.v)
  b = evaluate(L[2], P.v)
  c = evaluate(L[3], P.v)
  if iszero(a) && iszero(b) && iszero(c)
     error("This is a singular point of the curve")
  else
     V = gens(parent(C.eq))
     T = a*V[1] + b*V[2] + c*V[3]
     return ProjPlaneCurve(T)
  end
end

################################################################################
# Union of two projective plane curves (with multiplicity)

@doc Markdown.doc"""
    union(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem

Return the union of `C` and `D` (with multiplicity).
"""
function Base.union(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  F = C.eq*D.eq
  return ProjPlaneCurve(F)
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
    curve_components(C::ProjPlaneCurve)

Return a dictionary containing the irreducible components of `C` and their multiplicity.
"""
function curve_components(C::ProjPlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  return C.components
end

################################################################################
# Check irreducibility.
# TODO: change for a direct irreducibility check when available.

@doc Markdown.doc"""
    isirreducible(C::ProjPlaneCurve)

Return `true` if `C` is irreducible, and `false` otherwise.
"""
function Oscar.isirreducible(C::ProjPlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
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
    isreduced(C::ProjPlaneCurve)

Return `true` if `C` is reduced, and `false` otherwise.
"""
function Oscar.isreduced(C::ProjPlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  return all(isone, values(C.components))
end

################################################################################
# Compute squarefree equation.
# TODO: change for a direct squarefree computation when available.

@doc Markdown.doc"""
    reduction(C::ProjPlaneCurve)

Return the projective plane curve defined by the squarefree part of the equation of `C`.
"""
function reduction(C::ProjPlaneCurve)
  if isempty(C.components)
     f = factor(C.eq)
     C.components = Dict(ProjPlaneCurve(x) => f.fac[x] for x in keys(f.fac))
  end
  F = prod(D -> D.eq, keys(C.components))
  rC = ProjPlaneCurve(F)
  rC.components = Dict(ProjPlaneCurve(D.eq) => 1 for D in keys(C.components))
  return rC
end

################################################################################
# gives the common components of two projective plane curves

@doc Markdown.doc"""
    common_component(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem

Return the projective plane curve consisting of the common component of `C` and `D`, or an empty vector if they do not have a common component.
"""
function common_component(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq.f, D.eq.f)
  if isone(G)
     return Vector{ProjPlaneCurve}()
  else
     return ProjPlaneCurve(G)
  end
end

################################################################################
# dehomogenization with respect to the specified variable into the specified ring

@doc Markdown.doc"""
    dehomogenization([r::MPolyRing], F::Oscar.MPolyElem_dec, i::Int64)

Return the polynomial obtained by dehomogenization of `F` with respect to the `i`th variable of `parent(F)`. The new polynomial is in `r` if specified, or in a new ring with one variable less otherwise.
"""
function dehomogenization(r::MPolyRing{S}, F::Oscar.MPolyElem_dec{S}, i::Int64) where S <: FieldElem
  @assert ishomogenous(F)
  R = parent(F)
  if nvars(R) -1 != nvars(r)
     error("Incompatible number of variables")
  end
  V = gens(r)
  insert!(V, i, r(1))
  phi = hom(R.R, r, V)
  return phi(F)
end

################################################################################
# dehomogenization without specifying the ring with repect to the specified variable

function dehomogenization(F::Oscar.MPolyElem_dec, i::Int64)
  R = parent(F)
  A = String.(symbols(R))
  r = PolynomialRing(R.R.base_ring, deleteat!(A, i))
  dehomogenization(r[1], F, i)
end

################################################################################
# non decorated version

function dehomogenization(F::Oscar.MPolyElem, i::Int64)
  R = parent(F)
  A = grade(R)
  dehomogenization(A(F), i)
end

################################################################################
# non decorated version

function dehomogenization(r::MPolyRing{S}, F::Oscar.MPolyElem{S}, i::Int64) where S <: FieldElem
  R = parent(F)
  A = grade(R)
  dehomogenization(r, A(F), i)
end

###!#############################################################################
# convert array of lenght 2 to ProjSpcElem with 1 for last coordinate.
# Helping function

function Array_to_ProjSpcElem(PP::Oscar.Geometry.ProjSpc{S}, p::Array{S, 1}) where S <: FieldElem
  if dim(PP) != length(p)
     error("Not the right size")
  end
  m = push!(p, 1)
  return Oscar.Geometry.ProjSpcElem(PP, m)
end

################################################################################
# Compute the intersection of two curves. The output is a list: the first
# element the affine plane curve defined by the common components of the two
# curves (or is empty if no common component), the second element is the list
# of intersection points. Some might be contained in the common component too.

@doc Markdown.doc"""
    curve_intersect([PP::Oscar.Geometry.ProjSpc{S}], C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem

Return a list whose first element is the projective plane curve defined by the gcd of `C.eq` and `D.eq`, the second element is the list of the remaining intersection points when the common components are removed from `C` and `D` (the points are in `PP` if specified, or in a new projective space otherwise).
"""
function curve_intersect(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq.f, D.eq.f)
  R = parent(C.eq)
  CC = []
  if !isone(G)
     # We divide by the gcd to get curves without common components.
     F = R(div(C.eq.f, G))
     H = R(div(D.eq.f, G))
     push!(CC, ProjPlaneCurve(G))
  else
     F = C.eq
     H = D.eq
  end
  Pts = []
  r, (X, Y) = PolynomialRing(R.R.base_ring, ["X", "Y"])
  Fa = dehomogenization(r, F, 3)
  Ha = dehomogenization(r, H, 3)
  if !isconstant(Fa) && !isconstant(Ha)
     Ca = AffinePlaneCurveModule.AffinePlaneCurve(Fa)
     Da = AffinePlaneCurveModule.AffinePlaneCurve(Ha)
     L = AffinePlaneCurveModule.curve_intersect(Ca, Da)
     for p in L[2]
        push!(Pts, Array_to_ProjSpcElem(PP, p.coord))
     end
  end
  rr, (x,) = PolynomialRing(R.R.base_ring, ["x"])
  phi = hom(R.R, rr, [gens(rr)[1], rr(1), rr(0)])
  phiF = phi(F)
  phiH = phi(H)
  g = gcd(phiF, phiH)
  f = factor(g)
  ro = []
  for h in keys(f.fac)
     if total_degree(h) == 1
        f = h//lc(h)
        push!(ro, -f + gens(rr)[1])
     end
  end
  for y in ro
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [lc(y), R.R.base_ring(1), R.R.base_ring(0)]))
  end
  if iszero(evaluate(F, [1, 0, 0])) && iszero(evaluate(H, [1, 0, 0]))
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [R.R.base_ring(1), R.R.base_ring(0), R.R.base_ring(0)]))
  end
  return [CC,Pts]
end

################################################################################
# without specifying the projective space.

function curve_intersect(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
   R = parent(C.eq)
   PP = projective_space(R.R.base_ring, 2)
   curve_intersect(PP[1], C, D)
end

################################################################################
# Compute the singular locus of an affine plane curve by intersecting the
# derivatives. The points might also be contained in the components.

@doc Markdown.doc"""
    curve_singular_locus([PP::Oscar.Geometry.ProjSpc{S}], C::ProjPlaneCurve{S}) where S <: FieldElem

Return the reduced singular locus of `C` as a list whose first element is the projective plane curve consisting of the singular components of `C` (if any), and the second element is the list of the singular points of the reduction of `C` (the points are in `PP` if specified, or in a new projective space otherwise). The singular component might not contain any point over the considered field.
"""
function curve_singular_locus(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}) where S <: FieldElem
  D = reduction(C)
  Pts = Array{Oscar.Geometry.ProjSpcElem, 1}()
  CC = Array{ProjPlaneCurve, 1}()
  # The components with multiplicity > 1 are singular
  f = []
  for (h, c) in C.components
     if c != 1
        push!(f, h.eq)
     end
  end
  if !isempty(f)
     g = prod(f)
     push!(CC, ProjPlaneCurve(g))
  end
  J = jacobi_ideal(D)
  FX = gens(J)[1]
  FY = gens(J)[2]
  FZ = gens(J)[3]
  # Computation of the singular locus of the reduction of C.
  # With the fact that the equation is now squarefree, only points can appear.
  # We compute the singular points in the chart z=1, then the ones of the form
  # (x:1:0) and then if (1:0:0) is singular.
  Fa = dehomogenization(D.eq, 3)
  if !isconstant(Fa)
     Da = AffinePlaneCurve(Fa)
     L = AffinePlaneCurveModule.curve_singular_locus(Da)
     if !isempty(L[2])
        for p in L[2]
           push!(Pts, Array_to_ProjSpcElem(PP, p.coord))
        end
     end
  end
  R = parent(C.eq)
  rr, (x) = PolynomialRing(R.R.base_ring, ["x"])
  phi = hom(R.R, rr, [gens(rr)[1], rr(1), rr(0)])
  pF = phi(D.eq)
  pX = phi(FX)
  pY = phi(FY)
  pZ = phi(FZ)
  I = ideal([pF, pX, pY, pZ])
  g = groebner_basis(I, :lex, complete_reduction=true)
  f = factor(g[1])
  ro = []
  for h in keys(f.fac)
     if total_degree(h) == 1
        f = h//lc(h)
        push!(ro, -f + gens(rr)[1])
     end
  end
  for y in ro
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [lc(y), R.R.base_ring(1), R.R.base_ring(0)]))
  end
  if iszero(evaluate(D.eq, [1,0,0])) && iszero(evaluate(FX, [1,0,0])) && iszero(evaluate(FY, [1,0,0])) && iszero(evaluate(FZ, [1,0,0]))
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [R.R.base_ring(1), R.R.base_ring(0), R.R.base_ring(0)]))
  end
  return [CC, Pts]
end

################################################################################
# without specifying the projective space.

function curve_singular_locus(C::ProjPlaneCurve)
   R = parent(C.eq)
   PP = projective_space(R.R.base_ring, 2)
   curve_singular_locus(PP[1], C)
end

################################################################################
#

@doc Markdown.doc"""
    issmooth_curve(C::ProjPlaneCurve)

Return `true` if `C` has no singular point, and `false` otherwise.
"""
function issmooth_curve(C::ProjPlaneCurve)
   R = parent(C.eq)
   PP = projective_space(R.R.base_ring, 2)
   s = curve_singular_locus(PP[1], C)
   if isempty(s[1])
      return isempty(s[2])
   else
      error("The curve is not reduced.")
   end
end

################################################################################
end
using .ProjPlaneCurveModule
