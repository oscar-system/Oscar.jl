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

mutable struct ProjPlaneCurve{S}
  eq::Oscar.MPolyElem_dec{S}              # Equation of the curve (polynomial in three variables)
  degree::Int64                           # degree of the equation of the curve
  dimension::Int64                        # dimension of the curve (as a variety)
  components::Dict
  function ProjPlaneCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    r = new{S}()
    r.eq = eq
    if nvars(parent(r.eq)) != 3
       error("The defining equation must belong to a ring with three variables")
    elseif isconstant(r.eq) == true
       error("The defining equation must be non constant")
    else
       if ishomogenous(r.eq) == false
          error("The defining equation is not homogeneous")
       end
    end
    r.degree = -1                      # -1 when it is not computed yet
    r.dimension = 1                    # since C is a plane curve, the dimension is always 1
    r.components = Dict{ProjPlaneCurve{S}, Int64}()
    return r
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
> Returns the degree of the defining polynomial of C.
"""
function Oscar.degree(C::ProjPlaneCurve)
  if C.degree >= 0
     return C.degree
  else
     C.degree = total_degree(C.eq.f)
     return C.degree
  end
end

################################################################################
# Compute the Jacobian Ideal of C

@doc Markdown.doc"""
jacobi_ideal(C::ProjPlaneCurve)
> Returns the Jacobian ideal of the defining polynomial of C.
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
> Returns an error if P is not a point of C, false if P is a singular point of C, and true if P is a smooth point of C.
"""
function Oscar.issmooth(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  if dim(P.parent) != 2
     error("The point needs to be in a projective two dimensional space")
  elseif evaluate(C.eq, P.v) != 0
     error("The point is not on the curve defined by ", C.eq)
  else
     J = jacobi_ideal(C)
     L = gens(J)
     a = evaluate(L[1], P.v)
     b = evaluate(L[2], P.v)
     c = evaluate(L[3], P.v)
     if a == 0 && b == 0 && c == 0
        return false
     else
        return true
     end
  end
end

################################################################################
# If P is a smooth point of a projective plane curve C, compute the tangent of C
# at the point P.

@doc Markdown.doc"""
tangent(C::ProjPlaneCurve{S}, P::Point{S}) where S <: FieldElem
> Returns the tangent of C at P when P is a smooth point of C, and an error otherwise.
"""
function tangent(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  if dim(P.parent) != 2
     error("The point needs to be in a projective two dimensional space")
  elseif evaluate(C.eq, P.v) != 0
     error("The point is not on the curve")
  else
     J = jacobi_ideal(C)
     L = gens(J)
     a = evaluate(L[1], P.v)
     b = evaluate(L[2], P.v)
     c = evaluate(L[3], P.v)
     if a == 0 && b == 0 && c == 0
        error("This is a singular point of the curve")
     else
        V = gens(parent(C.eq))
        T = a*V[1] + b*V[2] + c*V[3]
        return ProjPlaneCurve(T)
     end
  end
end

################################################################################
# Union of two projective plane curves (with multiplicity)

@doc Markdown.doc"""
union(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
> Returns the union of C and D (with multiplicity).
"""
function Base.union(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  F = C.eq*D.eq
  return ProjPlaneCurve(F)
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
curve_components(C::ProjPlaneCurve)
> Returns a dictionary containing the irreducible components of C and their multiplicity.
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
> Returns true if C is irreducible, and false otherwise.
"""
function Oscar.isirreducible(C::ProjPlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  if length(C.components) != 1
     return false
  else
     return all(x -> x == 1, values(C.components))
  end
end

################################################################################
# Check reducedness by computing a factorization
# TODO: change for a direct squarefree check when available.

@doc Markdown.doc"""
isreduced(C::ProjPlaneCurve)
> Returns true if C is reduced, and false otherwise.
"""
function Oscar.isreduced(C::ProjPlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(ProjPlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  return all(x -> x == 1, values(C.components))
end

################################################################################
# Compute squarefree equation.
# TODO: change for a direct squarefree computation when available.

@doc Markdown.doc"""
reduction(C::ProjPlaneCurve)
> Returns the projective plane curve defined by the squarefree part of the equation of C.
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
> Returns the projective plane curve consisting of the common component of C and D, or an empty array if they do not have a common component.
"""
function common_component(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq.f, D.eq.f)
  if isone(G)
     return Array{ProjPlaneCurve, 1}()
  else
     return ProjPlaneCurve(G)
  end
end

################################################################################
# dehomogenization with respect to the specified variable into the specified ring

@doc Markdown.doc"""
dehomogenization(r::MPolyRing, F::Oscar.MPolyElem_dec, i::Int64)
> Returns the polynomial in r obtained by dehomogenization with respect to the ith variable of parent(F).
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

@doc Markdown.doc"""
dehomogenization(F::Oscar.MPolyElem_dec, i::Int64)
> Returns the polynomial in a new ring with one variable less obtained by dehomogenization with respect to the ith variable of parent(F).
"""
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
  q = Oscar.Geometry.ProjSpcElem(PP, m)
  return q
end

################################################################################
# Compute the intersection of two curves. The output is a list: the first
# element the affine plane curve defined by the common components of the two
# curves (or is empty if no common component), the second element is the list
# of intersection points. Some might be contained in the common component too.

@doc Markdown.doc"""
curve_intersect(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
> Returns a list whose first element is the projective plane curve defined by the gcd of F and G, the second element is the list of the remaining intersection points when the common components are removed from C and D (the points are in PP).
"""
function curve_intersect(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq.f, D.eq.f)
  R = parent(C.eq)
  CC = []
  if isone(G) == false
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
  if isconstant(Fa) == false && isconstant(Ha) == false
     Ca = AffinePlaneCurveModule.AffinePlaneCurve(Fa)
     Da = AffinePlaneCurveModule.AffinePlaneCurve(Ha)
     L = AffinePlaneCurveModule.curve_intersect(Ca, Da)
     for p in L[2]
        push!(Pts, Array_to_ProjSpcElem(PP, p.coord))
     end
  end
  rr, (x) = PolynomialRing(R.R.base_ring, ["x"])
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
  if evaluate(F, [1, 0, 0]) == 0 && evaluate(H, [1, 0, 0]) == 0
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [R.R.base_ring(1), R.R.base_ring(0), R.R.base_ring(0)]))
  end
  return [CC,Pts]
end

################################################################################
# without specifying the projective space.

@doc Markdown.doc"""
curve_intersect(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
> Returns a list whose first element is the projective plane curve defined by the gcd of F and G, the second element is the list of the remaining intersection points when the common components are removed from C and D (the points are in a new projective space).
"""
function curve_intersect(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
   R = parent(C.eq)
   PP = projective_space(R.R.base_ring, 2)
   curve_intersect(PP[1], C, D)
end

################################################################################
# Compute the singular locus of an affine plane curve by intersecting the
# derivatives. The points might also be contained in the components.

@doc Markdown.doc"""
curve_singular_locus(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}) where S <: FieldElem
> Returns a list whose first element is the projective plane curve consisting of the singular components of C (if any), and the second element is the list of the singular points of the reduction of C (the points are in PP). The singular component might not contain any point over the considered field.
"""
function curve_singular_locus(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}) where S <: FieldElem
  D = reduction(C)
  Pts = Array{Oscar.Geometry.ProjSpcElem, 1}()
  CC = Array{ProjPlaneCurve, 1}()
  # The components with multiplicity > 1 are singular
  f = []
  for h in keys(C.components)
     if C.components[h] != 1
        push!(f, h.eq^C.components[h])
     end
  end
  if isempty(f) == false
     g = prod(f)
     push!(CC, ProjPlaneCurve(g))
  end
  # Computation of the singular locus of the reduction of C (which consists of only isolated points)
  J = jacobi_ideal(D)
  FX = gens(J)[1]
  FY = gens(J)[2]
  FZ = gens(J)[3]
  if isconstant(FX) == true
     if FX != 0
        return [CC, Pts]
     else
        #FX = 0
        if isconstant(FY) == true
           if FY != 0
             return [CC, Pts]
           else
              #FX = FY = 0
              if isconstant(FZ) == true
                 if FZ != 0
                    return [CC, Pts]
                 else
                    # FX = FY = FZ = 0
                    push!(CC, C)
                    return [CC, Pts]
                 end
              else
                 # FX = FY = 0 and FZ non constant
                 return(CC, Pts)
              end
           end
        else
           # FX = 0 and FY non constant
           if isconstant(FZ) == true
              return [CC, Pts]
           else
              # FX = 0 and FY, FZ non constant
              CY = ProjPlaneCurve(FY)
              CZ = ProjPlaneCurve(FZ)
              L = curve_intersect(PP, CY, CZ)
              return [CC, L[2]]
           end
        end
     end
  else
     # FX non constant
     if isconstant(FY) == true
        if FY != 0
           return [CC, Pts]
        else
           # FX non constant, FY = 0
           if isconstant(FZ) == true
             return [CC, Pts]
           else
              # FX and FZ non constant, FY = 0
              CX = ProjPlaneCurve(FX)
              CZ = ProjPlaneCurve(FZ)
              L = curve_intersect(PP, CX, CZ)
              return [CC, L[2]]
           end
        end
     else
        # FX and FY non constant
        if isconstant(FZ) == true
           if FZ != 0
              return [CC, Pts]
           else
              # FX and FY non constant, FZ = 0
              CX = ProjPlaneCurve(FX)
              CY = ProjPlaneCurve(FY)
              L = curve_intersect(PP, CX, CY)
              return [CC, L[2]]
           end
        else
           # FX, FY, FZ non constant
           CX = ProjPlaneCurve(FX)
           CY = ProjPlaneCurve(FY)
           CZ = ProjPlaneCurve(FZ)
           L = curve_intersect(PP, CX, CY)
           for p in L[2]
              if evaluate(FZ, p.v) == 0
                 push!(Pts, p)
              end
           end
           if isempty(L[1])
              return [CC, Pts]
           else
              M = curve_intersect(PP, CZ, L[1][1])
              PP = [Pts; M[2]]
              return [CC, PP]
           end
        end
     end
  end
end

################################################################################
# without specifying the projective space.

@doc Markdown.doc"""
curve_singular_locus(C::ProjPlaneCurve)
> Returns a list whose first element is the projective plane curve consisting of the singular components of C (if any), and the second element is the list of the singular points of the reduction of C (the points are in a new projective space). The singular component might not contain any point over the considered field.
"""
function curve_singular_locus(C::ProjPlaneCurve)
   R = parent(C.eq)
   PP = projective_space(R.R.base_ring, 2)
   curve_singular_locus(PP[1], C)
end

################################################################################
#

@doc Markdown.doc"""
issmooth_curve(C::ProjPlaneCurve)
> Returns true if C has no singular point, and false otherwise.
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
