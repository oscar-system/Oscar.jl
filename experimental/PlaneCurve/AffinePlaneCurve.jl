
#module AffinePlaneCurveModule

include("Variety.jl")
using .VarietyModule

export issmooth, tangent, common_components, curve_intersect, intersect,
       curve_singular_locus, multiplicity, tangent_lines, singular_locus,
       intersection_multiplicity, aretransverse, issmooth_curve,
       arithmetic_genus, geometric_genus

################################################################################
# To check if a point is smooth: return true if the point is a smooth point on
# the curve C, and false  if it is singular, and an error if it is not on the
# curve.

@doc Markdown.doc"""
    issmooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem

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
# gives the common components of two affine plane curves

@doc Markdown.doc"""
    common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem(C::AffinePlaneCurve)

Return the affine plane curve consisting of the common component of `C` and `D`, or an empty vector if they do not have a common component.
"""
function common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  if isone(G)
     return Vector{AffinePlaneCurve}()
  else
     return [AffinePlaneCurve(G)]
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
  B = eliminate(I, [gen(R, 1)])
  # We factorize the polynomial we obtain (which is a polynomial where only
  # the second variable appears).
  Z = factor(gen(B, 1))
  Y = Vector{Oscar.MPolyElem}()
  L = Vector{S}[]
  # For the linear factors of g, we add the constant coefficient to the list Y.
  # The constant coefficient is stored as a mpoly to be used later in evaluate.
  for g in keys(Z.fac)
     if total_degree(g) == 1
        f = g//leading_coefficient(g)
        push!(Y, -f + gen(R, 2))
     end
  end
     # For each y, we compute the possible values of x by replacing the second
     # variable by the value y, and factorizing the resulting polynomial.
     for y in Y
        FF = evaluate(F, [gen(R, 1), y])
        HH = evaluate(H, [gen(R, 1), y])
        GG = gcd(FF, HH)
        ZZ = factor(GG)
        for g in keys(ZZ.fac)
           if total_degree(g) == 1
              f = g//leading_coefficient(g)
              # We use lc to convert the mpoly into element of S.
              push!(L, [leading_coefficient(-f + gen(R, 1)), leading_coefficient(y)])
           end
        end
     end
  # L contains the intersection points as arrays (or can be empty).
  Pts = [Point(p) for p in L]
  CC = Vector{AffinePlaneCurve}()
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
   comp = curve_components(C)
   D = reduction(C)
   Pts = Array{Point, 1}()
   CC = Array{AffinePlaneCurve, 1}()
   # The components with multiplicity > 1 are singular
   f = [h.eq for (h,c) in comp if c != 1]
   if !isempty(f)
      g = prod(f)
      push!(CC, AffinePlaneCurve(g))
   end
   # Computation of the singular locus of the reduction of C.
   J = jacobi_ideal(D)
   F = D.eq
   FX, FY = gens(J)
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
  A = grade(R)[1]
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
  A = grade(R)[1]
  HC = homogenous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  Gm = HC[M[1]]
  Z = factor(Gm.f)
  D = Dict{AffinePlaneCurve{S}, Int}()
  X = V[1] - P.coord[1]
  Y = V[2] - P.coord[2]
  for p in keys(Z.fac)
     if total_degree(p) == 1
        push!(D, AffinePlaneCurve(evaluate(p, [X, Y])) => Z.fac[p])
     end
  end
  return D
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

@doc Markdown.doc"""
    arithmetic_genus(C::AffinePlaneCurve)

Return the arithmetic genus of the projective closure of `C`.
"""
function arithmetic_genus(C::AffinePlaneCurve)
   F = defining_equation(C)
   G = homogenization(F, variable = "vrbl")
   D = ProjPlaneCurve(G)
   return arithmetic_genus(D)
end

################################################################################

@doc Markdown.doc"""
    geometric_genus(C::AffinePlaneCurve)

Return the geometric genus of the projective closure of `C`.
"""
function geometric_genus(C::AffinePlaneCurve)
   F = defining_equation(C)
   G = homogenization(F, variable = "vrbl")
   D = ProjPlaneCurve(G)
   return geometric_genus(D)
end

################################################################################
#end
#using .AffinePlaneCurveModule
