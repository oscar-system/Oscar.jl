module AffinePlaneCurveModule

using Oscar, Markdown
include("Variety.jl")
using .VarietyModule

import Base.:(==)

export AffinePlaneCurve, degree, jacobi_ideal, tangent, issmooth, union,
       Point, curve_components, isirreducible, isreduced, common_component,
       curve_intersect, intersect, curve_singular_locus, ideal_point,
       multiplicity, singular_locus, intersection_multiplicity,
       aretransverse, tangent_lines, issmooth_curve, reduction, hash

################################################################################
#
# We follow the lecture notes by Andreas Gathmann:
# https://www.mathematik.uni-kl.de/~gathmann/en/curves.php
#
################################################################################
#

mutable struct AffinePlaneCurve{S <: FieldElem}
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

Return the degree of the defining polynomial of `C`.
"""
function Oscar.degree(C::AffinePlaneCurve)
  if C.degree < 0
     C.degree = total_degree(C.eq)
  end
  return C.degree
end

################################################################################
# Compute the Jacobian Ideal of C

@doc Markdown.doc"""
    jacobi_ideal(C::AffinePlaneCurve)

Return the Jacobian ideal of the defining polynomial of `C`.
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

Throw an error if `P` is not a point of `C`, return `false` if `P` is a singular point of `C`, and `true` if `P` is a smooth point of `C`.
"""
function Oscar.issmooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  P.ambient_dim == 2 || error("The point needs to be in a two dimensional space")
  iszero(evaluate(C.eq, P.coord)) || error("The point is not on the curve defined by ", C.eq)
  J = jacobi_ideal(C)
  L = gens(J)
  a = evaluate(L[1], P.coord)
  b = evaluate(L[2], P.coord)
  return !(iszero(a) && iszero(b))
end

################################################################################
# If P is a smooth point of an affine plane curve C, compute the tangent of C
# at the point P.

@doc Markdown.doc"""
    tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem

Return the tangent of `C` at `P` when `P` is a smooth point of `C`, and throw an error otherwise.
"""
function tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  P.ambient_dim == 2 || error("The point needs to be in a two dimensional space")
  iszero(evaluate(C.eq, P.coord)) || error("The point is not on the curve")
  J = jacobi_ideal(C)
  L = gens(J)
  a = evaluate(L[1], P.coord)
  b = evaluate(L[2], P.coord)
  if iszero(a) && iszero(b)
     error("This is a singular point of the curve")
  else
     V = gens(parent(C.eq))
     T = a*(V[1] - P.coord[1]) + b*(V[2] - P.coord[2])
     return AffinePlaneCurve(T)
  end
end

################################################################################
# Union of two affine plane curves (with multiplicity)

@doc Markdown.doc"""
    union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem

Return the union of `C` and `D` (with multiplicity).
"""
function Base.union(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  F = C.eq*D.eq
  return AffinePlaneCurve(F)
end

################################################################################
# Components of the curve

@doc Markdown.doc"""
    curve_components(C::AffinePlaneCurve)

Return a dictionary containing the irreducible components of `C` and their multiplicity.
"""
function curve_components(C::AffinePlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  return C.components
end

################################################################################
# Check irreducibility.
# TODO: change for a direct irreducibility check when available.

@doc Markdown.doc"""
    isirreducible(C::AffinePlaneCurve)

Return `true` if `C` is irreducible, and `false` otherwise.
"""
function Oscar.isirreducible(C::AffinePlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
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
    isreduced(C::AffinePlaneCurve)

Return `true` if `C` is reduced, and `false` otherwise.
"""
function Oscar.isreduced(C::AffinePlaneCurve)
  if isempty(C.components)
     D = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => D.fac[x] for x in keys(D.fac))
  end
  return all(isone, values(C.components))
end

################################################################################
# gives the common components of two affine plane curves

@doc Markdown.doc"""
    common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem(C::AffinePlaneCurve)

Return the affine plane curve consisting of the common component of `C` and `D`, or an empty vector if they do not have a common component.
"""
function common_component(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  if isone(G)
     return Vector{AffinePlaneCurve}()
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

Return a list whose first element is the affine plane curve defined by the gcd of `C.eq` and `D.eq`, the second element is the list of the remaining intersection points when the common components are removed from `C` and `D`.
"""
function curve_intersect(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  R = parent(C.eq)
  if !isone(G)
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
  L = Array{S, 1}[]
  # For the linear factors of g, we add the constant coefficient to the list Y.
  # The constant coefficient is stored as a mpoly to be used later in evaluate.
  for g in keys(Z.fac)
     if total_degree(g) == 1
        f = g//lc(g)
        push!(Y, -f + gens(R)[2])
     end
  end
  if !isempty(Y)
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
  if !isempty(L)
     for p in L
        push!(Pts, Point(p))
     end
  end
  if !isone(G)
     push!(CC, AffinePlaneCurve(G))
  end
  return [CC, Pts]
end

################################################################################
# Intersection in the sense of varieties.
# Might change depending on future changes on VarietyModule.

@doc Markdown.doc"""
    intersect( C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem

Return the variety defined by the intersection of `C` and `D`.
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
    curve_singular_locus(C::AffinePlaneCurve)

Return the reduced singular locus of `C` as a list whose first element is the affine plane curve consisting of the singular components of `C` (if any), and the second element is the list of the isolated singular points (which may be contained in the singular component). The singular component might not contain any point over the considered field.
"""
function curve_singular_locus(C::AffinePlaneCurve)
   D = reduction(C)
   Pts = Array{Point, 1}()
   CC = Array{AffinePlaneCurve, 1}()
   # The components with multiplicity > 1 are singular
   f = []
   for (h, c) in C.components
      if c != 1
         push!(f, h.eq)
      end
   end
   if !isempty(f)
      g = prod(f)
      push!(CC, AffinePlaneCurve(g))
   end
   # Computation of the singular locus of the reduction of C.
   J = jacobi_ideal(D)
   F = D.eq
   FX = gens(J)[1]
   FY = gens(J)[2]
   # The case FX = FY = 0 cannot occur in the reduced case.
   # With the fact that the equation is now squarefree, only points can appear.
   if iszero(FX) && !isconstant(FY)
      CY = AffinePlaneCurve(FY)
      L = curve_intersect(D, CY)
      append!(Pts, L[2])
   elseif iszero(FY) && !isconstant(FX)
      CX = AffinePlaneCurve(FX)
      L = curve_intersect(D, CX)
      append!(Pts, L[2])
   elseif !isconstant(FX) && !isconstant(FY)
      CX = AffinePlaneCurve(FX)
      CY = AffinePlaneCurve(FY)
      L = curve_intersect(CX, CY)
      if !isempty(L[1])
         M = curve_intersect(L[1][1], D)
         append!(Pts, M[2])
      end
      if !isempty(L[2])
         for p in L[2]
            if iszero(evaluate(F, p.coord))
               push!(Pts, p)
            end
         end
      end
   end
   return [CC, Pts]
end

################################################################################
# Compute squarefree equation.
# TODO: change for a direct squarefree computation when available.

@doc Markdown.doc"""
    reduction(C::AffinePlaneCurve)

Return the affine plane curve defined by the squarefree part of the equation of `C`.
"""
function reduction(C::AffinePlaneCurve)
  if isempty(C.components)
     f = factor(C.eq)
     C.components = Dict(AffinePlaneCurve(x) => f.fac[x] for x in keys(f.fac))
  end
  F = prod(D -> D.eq, keys(C.components))
  rC = AffinePlaneCurve(F)
  rC.components = Dict(AffinePlaneCurve(D.eq) => 1 for D in keys(C.components))
  return rC
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
# Helping function:
# change of variable to send a point at the origin.

function curve_map_point_origin(C::AffinePlaneCurve, P::Point)
  F = C.eq
  R = parent(C.eq)
  V = gens(R)
  G = evaluate(F, [V[1] + P.coord[1], V[2] + P.coord[2]])
  return AffinePlaneCurve(G)
end

################################################################################
# Helping function:
# used to order elements to find the multiplicity.

function _sort_helper_multiplicity(a::GrpAbFinGenElem)
  return a.coeff[1, 1]
end

################################################################################
# compute the multiplicity of the affine plane curve C at the point P

@doc Markdown.doc"""
    multiplicity(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem

Return the multiplicity of `C` at `P`.
"""
function multiplicity(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  P.ambient_dim == 2 || error("The point needs to be in a two dimensional space")
  D = curve_map_point_origin(C, P)
  G = D.eq
  R = parent(G)
  A = grade(R)
  HC = homogenous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  return M[1].coeff[1, 1]
end

################################################################################
# compute the set of tangent lines of the affine plane curve C at the point P
# (linear factors of the homogeneous part of lower degree of the equation).

@doc Markdown.doc"""
    tangent_lines(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem

Return the tangent lines at `P` to `C` with their multiplicity.
"""
function tangent_lines(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  P.ambient_dim == 2 || error("The point needs to be in a two dimensional space")
  D = curve_map_point_origin(C, P)
  G = D.eq
  R = parent(G)
  V = gens(R)
  A = grade(R)
  HC = homogenous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  Gm = HC[M[1]]
  Z = factor(Gm)
  D = Dict{AffinePlaneCurve{S}, Int64}()
  X = V[1] - P.coord[1]
  Y = V[2] - P.coord[2]
  for p in keys(Z.fac)
     if total_degree(p.f) == 1
        push!(D, AffinePlaneCurve(evaluate(p.f, [X, Y])) => Z.fac[p])
     end
  end
  return D
end

################################################################################
# check for equality of curves (the equations are equal up to multiplication by
# a non zero constant).

function ==(C::AffinePlaneCurve, D::AffinePlaneCurve)
  F = C.eq
  G = D.eq
  return degree(C) == degree(D) && F*(lc(G)//lc(F)) == G
end

################################################################################
# Singular locus in the sense of varieties.
# Might change depending on future changes on VarietyModule.

@doc Markdown.doc"""
    singular_locus(C::AffinePlaneCurve)

Return the singular locus of `C` as a variety.
"""
function singular_locus(C::AffinePlaneCurve)
  J = jacobi_ideal(C)
  R = parent(C.eq)
  S = ideal(R, [gens(J); C.eq])
  return Variety(S)
end

################################################################################
# Compute the intersection multiplicity of the two affine plane curves at the
# given point.

@doc Markdown.doc"""
    intersection_multiplicity(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem

Return the intersection multiplicity of `C` and `D` at `P`.
"""
function intersection_multiplicity(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
  P.ambient_dim == 2 || error("The point needs to be in a two dimensional space")
  R = parent(C.eq)
  m = ideal_point(R, P)
  r = Localization(R, m)
  I = ideal(r, [C.eq, D.eq])
  groebner_basis(I)
  return Singular.vdim(I.gb.S)
end

################################################################################
# Check if the two curves intersect transversally at the given point.

@doc Markdown.doc"""
    aretransverse(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S<:FieldElem

Return `true` if `C` and `D` intersect transversally at `P` and `false` otherwise.
"""
function aretransverse(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S<:FieldElem
  return issmooth(C, P) && issmooth(D, P) && intersection_multiplicity(C, D, P) == 1
end

################################################################################
# Check if a reduced curve is smooth.

@doc Markdown.doc"""
    issmooth_curve(C::AffinePlaneCurve)

Return `true` if `C` has no singular point, and `false` otherwise.
"""
function issmooth_curve(C::AffinePlaneCurve)
  S = curve_singular_locus(C)
  if isempty(S[1])
     return isempty(S[2])
  else
     error("The curve is not reduced.")
  end
end

################################################################################
# hash function

function Base.hash(C::AffinePlaneCurve, h::UInt)
  F = 1//lc(C.eq)*C.eq
  return hash(F, h)
end

################################################################################

end
using .AffinePlaneCurveModule
