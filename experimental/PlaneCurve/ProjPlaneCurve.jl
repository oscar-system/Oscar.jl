#module ProjPlaneCurveModule

export is_smooth, tangent, common_components, curve_intersect,
       curve_singular_locus, is_smooth_curve, multiplicity,
       tangent_lines, intersection_multiplicity, aretransverse,
       arithmetic_genus, geometric_genus


################################################################################
# Functions for ProjPlaneCurves
################################################################################
# To check if a point is smooth: return true if the point is a smooth point on
# the curve C, and false  if it is singular, and an error if it is not on the
# curve.

@doc Markdown.doc"""
    is_smooth(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Throw an error if `P` is not a point of `C`, return `false` if `P` is a singular point of `C`, and `true` if `P` is a smooth point of `C`.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(x^2*(x+y)*(y^3-x^2*z))
Projective plane curve defined by -x^5*z - x^4*y*z + x^3*y^3 + x^2*y^4

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
(0 : 0 : 1)

julia> Oscar.is_smooth(C, P)
false
```
"""
function Oscar.is_smooth(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
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
    tangent(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the tangent of `C` at `P` when `P` is a smooth point of `C`, and throw an error otherwise.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y","z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> C = Oscar.ProjPlaneCurve(x^2*(x+y)*(y^3-x^2*z))
Projective plane curve defined by -x^5*z - x^4*y*z + x^3*y^3 + x^2*y^4

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
(2 : -2 : 1)

julia> Oscar.tangent(C, P)
Projective plane curve defined by -48*x - 48*y
```
"""
function tangent(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
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
    common_components(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}) where S <: FieldElem

Return the projective plane curve consisting of the common component of `C` and `D`, or an empty vector if they do not have a common component.
"""
function common_components(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}) where S <: FieldElem
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

function Array_to_ProjSpcElem(PP::Oscar.Geometry.ProjSpc{S}, p::Vector{S}) where S <: FieldElem
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
    curve_intersect([PP::Oscar.Geometry.ProjSpc{S}], C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}) where S <: FieldElem

Return a list whose first element is the projective plane curve defined by the gcd of `C.eq` and `D.eq`, the second element is the list of the remaining intersection points when the common components are removed from `C` and `D` (the points are in `PP` if specified, or in a new projective space otherwise).

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y","z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> C = Oscar.ProjPlaneCurve(T(x+y+z))
Projective plane curve defined by x + y + z

julia> D = Oscar.ProjPlaneCurve(T(z))
Projective plane curve defined by z

julia> Oscar.curve_intersect(PP[1], C, D)
2-element Vector{Vector{Any}}:
 []
 [(-1 : 1 : 0)]
```
"""
function curve_intersect(PP::Oscar.Geometry.ProjSpc{S}, C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}) where S <: FieldElem
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
  Fa = dehomogenization(F, r, 3)
  Ha = dehomogenization(H, r, 3)
  if !is_constant(Fa) && !is_constant(Ha)
     Ca = AffinePlaneCurve(Fa)
     Da = AffinePlaneCurve(Ha)
     L = curve_intersect(Ca, Da)
     for p in L[2]
        push!(Pts, Array_to_ProjSpcElem(PP, p.coord))
     end
  end
  rr, (x,) = PolynomialRing(R.R.base_ring, ["x"])
  phi = hom(R.R, rr, [gen(rr, 1), rr(1), rr(0)])
  phiF = phi(F.f)
  phiH = phi(H.f)
  g = gcd(phiF, phiH)
  f = factor(g)
  ro = []
  for h in keys(f.fac)
     if total_degree(h) == 1
        f = h//leading_coefficient(h)
        push!(ro, -f + gen(rr, 1))
     end
  end
  for y in ro
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [leading_coefficient(y), R.R.base_ring(1), R.R.base_ring(0)]))
  end
  if iszero(evaluate(F, [1, 0, 0])) && iszero(evaluate(H, [1, 0, 0]))
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [R.R.base_ring(1), R.R.base_ring(0), R.R.base_ring(0)]))
  end
  return [CC,Pts]
end

################################################################################
# without specifying the projective space.

function curve_intersect(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}) where S <: FieldElem
   R = parent(C.eq)
   PP = proj_space(R.R.base_ring, 2)
   curve_intersect(PP[1], C, D)
end

################################################################################
# Compute the singular locus of an affine plane curve by intersecting the
# derivatives. The points might also be contained in the components.

@doc Markdown.doc"""
    curve_singular_locus([PP::Oscar.Geometry.ProjSpc{S}], C::ProjectivePlaneCurve{S}) where S <: FieldElem

Return the reduced singular locus of `C` as a list whose first element is the projective plane curve consisting of the singular components of `C` (if any), and the second element is the list of the singular points of the reduction of `C` (the points are in `PP` if specified, or in a new projective space otherwise). The singular component might not contain any point over the considered field.
"""
function curve_singular_locus(PP::Oscar.Geometry.ProjSpc{S}, C::ProjectivePlaneCurve{S}) where S <: FieldElem
  D = reduction(C)
  Pts = Vector{Oscar.Geometry.ProjSpcElem}()
  CC = Vector{ProjPlaneCurve}()
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
  if !is_constant(Fa)
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
  pF = phi(D.eq.f)
  pX = phi(FX.f)
  pY = phi(FY.f)
  pZ = phi(FZ.f)
  I = ideal([pF, pX, pY, pZ])
  g = collect(Oscar.groebner_assure(I, lex(gens(rr)), true))
  f = factor(g[1])
  ro = []
  for h in keys(f.fac)
     if total_degree(h) == 1
        f = h//leading_coefficient(h)
        push!(ro, -f + gen(rr, 1))
     end
  end
  for y in ro
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [leading_coefficient(y), R.R.base_ring(1), R.R.base_ring(0)]))
  end
  if iszero(evaluate(D.eq, [1,0,0])) && iszero(evaluate(FX, [1,0,0])) && iszero(evaluate(FY, [1,0,0])) && iszero(evaluate(FZ, [1,0,0]))
     push!(Pts, Oscar.Geometry.ProjSpcElem(PP, [R.R.base_ring(1), R.R.base_ring(0), R.R.base_ring(0)]))
  end
  return [CC, Pts]
end

################################################################################
# without specifying the projective space.

function curve_singular_locus(C::ProjectivePlaneCurve)
   R = parent(C.eq)
   PP = proj_space(R.R.base_ring, 2)
   curve_singular_locus(PP[1], C)
end

################################################################################
#

@doc Markdown.doc"""
    is_smooth_curve(C::ProjectivePlaneCurve)

Return `true` if `C` has no singular point, and `false` otherwise.
"""
function is_smooth_curve(C::ProjectivePlaneCurve)
   F = defining_equation(C)
   R = parent(F)
   J = jacobi_ideal(F)
   if dim(J) > 0
      return false
   else
      return true
   end
end

################################################################################
# multiplicity

@doc Markdown.doc"""
     multiplicity(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the multiplicity of `C` at `P`.
"""
function multiplicity(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
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
      tangent_lines(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the tangent lines at `P` to `C` with their multiplicity.
 """
 function tangent_lines(C::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
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

function _dehom_curves_r(r::MPolyRing, C::ProjectivePlaneCurve, D::ProjectivePlaneCurve, i::Int)
  F = dehomogenization(C.eq, r, i)
  G = dehomogenization(D.eq, r, i)
  return [AffinePlaneCurve(F), AffinePlaneCurve(G)]
end

################################################################################

@doc Markdown.doc"""
     intersection_multiplicity(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the intersection multiplicity of `C` and `D` at `P`.
"""
function intersection_multiplicity(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
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
     aretransverse(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S<:FieldElem

Return `true` if `C` and `D` intersect transversally at `P` and `false` otherwise.
"""
function aretransverse(C::ProjectivePlaneCurve{S}, D::ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
  return is_smooth(C, P) && is_smooth(D, P) && intersection_multiplicity(C, D, P) == 1
end

################################################################################

@doc Markdown.doc"""
    arithmetic_genus(C::ProjectivePlaneCurve)

Return the arithmetic genus of `C`.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> C = Oscar.ProjPlaneCurve(T(y^2 * z - x^3 - x * z^2))
Projective plane curve defined by -x^3 - x*z^2 + y^2*z

julia> Oscar.arithmetic_genus(C)
1
```
"""
function arithmetic_genus(C::ProjectivePlaneCurve)
  F = C.eq
  S = parent(F)
  I = ideal(S, [F])
  A = quo(S, I)
  H = hilbert_polynomial(A[1])
  return -ZZ(coeff(H, 0)) + 1
end

################################################################################

@doc Markdown.doc"""
    geometric_genus(C::ProjectivePlaneCurve{T}) where T <: FieldElem

Return the geometric genus of `C`.

# Examples
```jldoctest
julia> R, (x,y,z) = GradedPolynomialRing(QQ, ["x", "y", "z"]);

julia> C = ProjPlaneCurve(z*x^2-y^3)
Projective plane curve defined by x^2*z - y^3

julia> geometric_genus(C)
0
```
"""
function geometric_genus(C::ProjectivePlaneCurve{S}) where S <: FieldElem
   F = defining_equation(C)
   R = parent(F)
   I = ideal(R, [F])
   if S == fmpq
      singular_assure(I)
      return ZZ(Singular.LibNormal.genus(I.gens.S)::Int)
   else
      A = quo(R, I)
      L = normalization(A[1])
      m = length(L)
      pa = zero(ZZ)
      for i in 1:m
         J = L[i][1].I
         T, _ = grade(parent(J[1]))
         V = T.(gens(J))
         JJ = ideal(T, V)
         B = quo(T, JJ)
         H = hilbert_polynomial(B[1])
         pa = pa - ZZ(coeff(H, 0))
      end
      return pa + 1
   end
end

################################################################################
#end
#using .ProjPlaneCurveModule
