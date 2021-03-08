#module ProjPlaneCurveModule

export dehomogenization, homogenization, issmooth, tangent,
       common_components, curve_intersect, curve_singular_locus, issmooth_curve,
       multiplicity, tangent_lines, intersection_multiplicity, aretransverse

################################################################################
# Homogenization and dehomogenization (not specific to curves).
################################################################################
# dehomogenization with respect to the specified variable into the specified
# ring

@doc Markdown.doc"""
    dehomogenization([r::MPolyRing], F::Oscar.MPolyElem_dec, i::Int)

Return the polynomial obtained by dehomogenization of `F` with respect to the `i`th variable of `parent(F)`. The new polynomial is in `r` if specified, or in a new ring with one variable less otherwise.
"""
function dehomogenization(r::MPolyRing{S}, F::Oscar.MPolyElem_dec{S}, i::Int) where S <: FieldElem
  @assert ishomogenous(F)
  R = parent(F)
  nvars(R) -1 == nvars(r) || error("Incompatible number of variables")
  V = gens(r)
  insert!(V, i, r(1))
  phi = hom(R.R, r, V)
  return phi(F)
end

################################################################################
# dehomogenization without specifying the ring with repect to the specified variable

function dehomogenization(F::Oscar.MPolyElem_dec, i::Int)
  R = parent(F)
  A = String.(symbols(R))
  r = PolynomialRing(R.R.base_ring, deleteat!(A, i))
  dehomogenization(r[1], F, i)
end

################################################################################
# non decorated version

function dehomogenization(F::Oscar.MPolyElem, i::Int)
  R = parent(F)
  A = grade(R)
  dehomogenization(A(F), i)
end

################################################################################
# non decorated version

function dehomogenization(r::MPolyRing{S}, F::Oscar.MPolyElem{S}, i::Int) where S <: FieldElem
  R = parent(F)
  A = grade(R)
  dehomogenization(r, A(F), i)
end

################################################################################
# homogenization

 @doc Markdown.doc"""
     homogenization(R::MPolyRing{S}, F::MPolyElem{S}, i::Int) where S <: FieldElem

Return the homogenization of `F` in `R` of a polynomial using the `i`th variable of `R`.
 """
 function homogenization(R::MPolyRing{S}, F::MPolyElem{S}, i::Int) where S <: FieldElem
   r = parent(F)
   V = gens(R)
   W = [V[j]//V[i] for j=1:nvars(R)]
   deleteat!(V, i)
   phi = hom(r, R, V)
   G = phi(F)
   d = total_degree(F)
   P = R[i]^d*evaluate(G, W)
   return numerator(P)
 end

################################################################################

@doc Markdown.doc"""
      homogenization(F::MPolyElem, x0::String="x0")

Return the homogenization of `F` in a ring with one additional variable (named `x0` if not specified).
 """
 function homogenization(F::MPolyElem, x0::String="x0")
   r = parent(F)
   A = String.(symbols(r))
   push!(A, x0)
   R = PolynomialRing(r.base_ring, A)
   homogenization(R[1], F, nvars(r)+1)
 end

################################################################################
# Functions for ProjPlaneCurves
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
    tangent(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

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
# gives the common components of two projective plane curves

@doc Markdown.doc"""
    common_components(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem

Return the projective plane curve consisting of the common component of `C` and `D`, or an empty vector if they do not have a common component.
"""
function common_components(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}) where S <: FieldElem
  G = gcd(C.eq, D.eq)
  if isone(G)
     return Vector{ProjPlaneCurve}()
  else
     return [ProjPlaneCurve(G)]
  end
end


################################################################################
# convert array of lenght 2 to ProjSpcElem with 1 for last coordinate.
# Helping function

function Array_to_ProjSpcElem(PP::Oscar.Geometry.ProjSpc{S}, p::Array{S, 1}) where S <: FieldElem
  dim(PP) == length(p) || error("Not the right size")
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
  G = gcd(C.eq, D.eq)
  R = parent(C.eq)
  CC = []
  if !isone(G)
     # We divide by the gcd to get curves without common components.
     F = div(C.eq, G)
     H = div(D.eq, G)
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
     Ca = AffinePlaneCurve(Fa)
     Da = AffinePlaneCurve(Ha)
     L = curve_intersect(Ca, Da)
     for p in L[2]
        push!(Pts, Array_to_ProjSpcElem(PP, p.coord))
     end
  end
  rr, (x,) = PolynomialRing(R.R.base_ring, ["x"])
  phi = hom(R.R, rr, [gen(rr, 1), rr(1), rr(0)])
  phiF = phi(F)
  phiH = phi(H)
  g = gcd(phiF, phiH)
  f = factor(g)
  ro = []
  for h in keys(f.fac)
     if total_degree(h) == 1
        f = h//lc(h)
        push!(ro, -f + gen(rr, 1))
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
  FX, FY, FZ = gens(J)
  # Computation of the singular locus of the reduction of C.
  # With the fact that the equation is now squarefree, only points can appear.
  # We compute the singular points in the chart z=1, then the ones of the form
  # (x:1:0) and then if (1:0:0) is singular.
  Fa = dehomogenization(D.eq, 3)
  if !isconstant(Fa)
     Da = AffinePlaneCurve(Fa)
     L = curve_singular_locus(Da)
     if !isempty(L[2])
        for p in L[2]
           push!(Pts, Array_to_ProjSpcElem(PP, p.coord))
        end
     end
  end
  R = parent(C.eq)
  rr, (x) = PolynomialRing(R.R.base_ring, ["x"])
  phi = hom(R.R, rr, [gen(rr, 1), rr(1), rr(0)])
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
        push!(ro, -f + gen(rr, 1))
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
# multiplicity

@doc Markdown.doc"""
     multiplicity(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Returns the multiplicity of `C` at `P`.
"""
function multiplicity(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  if P.v[3] != 0
     Fa = dehomogenization(C.eq, 3)
     Ca = AffinePlaneCurve(Fa)
     Q = Point([P.v[1]//P.v[3], P.v[2]//P.v[3]])
  elseif P.v[2] != 0
     Fa = dehomogenization(C.eq, 2)
     Ca = AffinePlaneCurve(Fa)
     Q = Point([P.v[1]//P.v[2], P.v[3]//P.v[2]])
  else
     Fa = dehomogenization(C.eq, 1)
     Ca = AffinePlaneCurve(Fa)
     Q = Point([P.v[2]//P.v[1], P.v[3]//P.v[1]])
  end
  return multiplicity(Ca, Q)
end

################################################################################
# homogeneization for lines

function help_homogene_line(R::MPolyRing, r::MPolyRing, F::MPolyElem, i::Int)
  total_degree(F) == 1 || error("This is not a degree one polynomial")
  V = gens(R)
  W = gens(r)
  v = V[i]
  deleteat!(V, i)
  phi = hom(r, R, V)
  G = phi(F)
  G = G - evaluate(G, [0,0,0])*(1-v)
  return G
end

################################################################################
# tangent lines

 @doc Markdown.doc"""
      tangent_lines(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Returns the tangent lines at `P` to `C` with their multiplicity.
 """
 function tangent_lines(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
   evaluate(C.eq, P.v) == 0 || error("The point is not on the curve.")
   R = parent(C.eq)
   if P.v[3] != 0
      Fa = dehomogenization(C.eq, 3)
      Ca = AffinePlaneCurve(Fa)
      Q = Point([P.v[1]//P.v[3], P.v[2]//P.v[3]])
      i = 3
   elseif P.v[2] != 0
      Fa = dehomogenization(C.eq, 2)
      Ca = AffinePlaneCurve(Fa)
      Q = Point([P.v[1]//P.v[2], P.v[3]//P.v[2]])
      i = 2
   else
      Fa = dehomogenization(C.eq, 1)
      Ca = AffinePlaneCurve(Fa)
      Q = Point([P.v[2]//P.v[1], P.v[3]//P.v[1]])
      i = 1
   end
   L = tangent_lines(Ca, Q)
   D = Dict{ProjPlaneCurve{S}, Int}()
   if isempty(L) == false
      r = parent(Ca.eq)
      D = Dict(ProjPlaneCurve(help_homogene_line(R, r, x.eq, i)) => L[x] for x in keys(L))
      #D = Dict(help_homogene_line(R, r, x.eq, i) => L[x] for x in keys(L))
   end
   return D
 end

################################################################################
# helping function for intersection_multiplicity

function _dehom_curves_r(r::MPolyRing, C::ProjPlaneCurve, D::ProjPlaneCurve, i::Int)
  F = dehomogenization(r, C.eq, i)
  G = dehomogenization(r, D.eq, i)
  return [AffinePlaneCurve(F), AffinePlaneCurve(G)]
end

################################################################################

@doc Markdown.doc"""
     intersection_multiplicity(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Returns the intersection multiplicity of `C` and `D` at `P`.
"""
function intersection_multiplicity(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
   dim(P.parent) == 2 || error("The point needs to be in a projective two dimensional space")
   R = parent(C.eq)
   r, (X, Y) = PolynomialRing(R.R.base_ring, ["X", "Y"])
   if P.v[3] != 0
      V = _dehom_curves_r(r, C, D, 3)
      Q = Point([P.v[1]//P.v[3], P.v[2]//P.v[3]])
   elseif P.v[2] != 0
      V = _dehom_curves_r(r, C, D, 2)
      Q = Point([P.v[1]//P.v[2], P.v[3]//P.v[2]])
   else
      V = _dehom_curves_r(r, C, D, 1)
      Q = Point([P.v[2]//P.v[1], P.v[3]//P.v[1]])
   end
   return intersection_multiplicity(V[1], V[2], Q)
end

################################################################################
#

@doc Markdown.doc"""
     aretransverse(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S<:FieldElem

Returns `true` if `C` and `D` intersect transversally at `P` and `false` otherwise.
"""
function aretransverse(C::ProjPlaneCurve{S}, D::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  return issmooth(C, P) && issmooth(D, P) && intersection_multiplicity(C, D, P) == 1
end

################################################################################
#end
#using .ProjPlaneCurveModule
