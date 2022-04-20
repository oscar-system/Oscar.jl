export ProjEllipticCurve, discriminant, issmooth, j_invariant,
       Point_EllCurve, curve, weierstrass_form, toweierstrass,
       iselliptic, list_rand

################################################################################
# Helping functions
################################################################################

function isweierstrass_form(eq::Oscar.MPolyElem{S}) where S <: RingElem
   R = parent(eq)
   x = gen(R, 1)
   y = gen(R, 2)
   z = gen(R, 3)
   A = [x*y*z, y*z^2, x^2*z, x*z^2, z^3]
   f = eq - y^2*z + x^3
   for t in A
      f = f - coeff(eq, t)*t
   end
   return iszero(f)
end

################################################################################

function shortformtest(eq::Oscar.MPolyElem{S}) where S <: RingElem
   R = parent(eq)
   x = gen(R, 1)
   y = gen(R, 2)
   z = gen(R, 3)
   A = [coeff(eq, x*y*z), -coeff(eq, x^2*z), coeff(eq, y*z^2), -coeff(eq, x*z^2), -coeff(eq, z^3)]
   iszero(A[1]) && iszero(A[2]) && iszero(A[3]) && return (true, [A[4], A[5]])
   return (false, A)
end

################################################################################

function _iselliptic(F::Oscar.MPolyElem_dec)
   C = ProjPlaneCurve(F)
   return iselliptic(C)
end

function iselliptic(C::ProjPlaneCurve)
   return degree(C) == 3 && issmooth_curve(C)
end

################################################################################

function isinflection(F::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
   J = jacobi_matrix([jacobi_matrix(F)[i] for i in 1:3])
   H = det(J)
   return iszero(evaluate(F, P.v)) && iszero(evaluate(H, P.v))
end

################################################################################
# Change of coordinates to get the Weierstrass form.
################################################################################

function _point_line(L::ProjPlaneCurve, Q::Oscar.Geometry.ProjSpcElem)
   degree(L) == 1 || error("Not a line")
   Q in L || error("the point is not on the line")
   F = L.eq
   R = parent(F)
   x = gen(R, 1)
   y = gen(R, 2)
   z = gen(R, 3)
   PP = Q.parent
   K = R.R.base_ring
   if coeff(F, gen(R, 3)) != 0
      if Q.v[2] != 0
         return Oscar.Geometry.ProjSpcElem(PP, [coeff(F, z), K(0), - coeff(F, x)])
      else
         return Oscar.Geometry.ProjSpcElem(PP, [K(0), coeff(F, z), - coeff(F, y)])
      end
   else
      if Q.v[1] != 0
         return Oscar.Geometry.ProjSpcElem(PP, [K(0), K(0), K(1)])
      else
         return Oscar.Geometry.ProjSpcElem(PP, [coeff(F, y), - coeff(F, x), K(0)])
      end
   end
end

################################################################################

function _completion(R::MPolyRing, P::Oscar.Geometry.ProjSpcElem, Q::Oscar.Geometry.ProjSpcElem)
   A = hcat(Q.v, P.v)
   if det(Oscar.matrix(R, 3, 3, hcat(A, [1, 0, 0]))) != 0
      B = Oscar.matrix(R, 3, 3, hcat(A, [1, 0, 0]))
   elseif det(Oscar.matrix(R, 3, 3, hcat(A, [0, 1, 0]))) != 0
      B = Oscar.matrix(R, 3, 3, hcat(A, [0, 1, 0]))
   else
      B = Oscar.matrix(R, 3, 3, hcat(A, [0, 0, 1]))
   end
   return [B, inv(B)]
end

################################################################################

function _change_coord(F::Oscar.MPolyElem_dec, P::Oscar.Geometry.ProjSpcElem)
   evaluate(F, P) == 0 || error("The point is not on the curve")
   D = ProjPlaneCurve(F)
   L = tangent(D, P)
   isinflection(F, P) || error("The point is not an inflection point")
   R = parent(defining_equation(D))
   S = parent(F)
   Q = _point_line(L, P)
   A = _completion(R, P, Q)
   V = Oscar.matrix(R, 3, 1, gens(R))
   B1 = A[1]*V
   C1 = [S(B1[i]) for i in 1:3]
   B2 = A[2]*V
   C2 = [S(B2[i]) for i in 1:3]
   return [hom(S, S, C1), hom(S, S, C2)]
end

################################################################################

function _discr(eq::Oscar.MPolyElem_dec)
   R = parent(eq)
   x = gen(R, 1)
   z = gen(R, 3)
   a = -coeff(eq, x*z^2)
   b = -coeff(eq, z^3)
   d = 4*a^3 + 27*b^2
   return d
end

################################################################################
################################################################################
# Definition EllipticCurve
################################################################################

################################################################################
# Definition EllipticCurve
################################################################################
@doc Markdown.doc"""
    ProjEllipticCurve{S}(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    ProjEllipticCurve(eq::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where {S <: FieldElem}
    ProjEllipticCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: Nemo.fmpz_mod}

Return the Projective Elliptic Curve defined by the equation `eq`, with `P` as
infinity point. If no point is specified it is expected that `eq` is in
Weierstrass form, and the infinity point is `(0:1:0)`.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> F = T(-x^3 - 3*x^2*y - 3*x*y^2 - x*z^2 - y^3 + y^2*z - y*z^2 - 4*z^3)
-x^3 - 3*x^2*y - 3*x*y^2 - x*z^2 - y^3 + y^2*z - y*z^2 - 4*z^3

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
(-1 : 1 : 0)

julia> E1 = Oscar.ProjEllipticCurve(F, P)
Projective elliptic curve defined by -x^3 - 3*x^2*y - 3*x*y^2 - x*z^2 - y^3 + y^2*z - y*z^2 - 4*z^3

julia> E2 = Oscar.ProjEllipticCurve(T(y^2*z - x^3 - x*z^2))
Projective elliptic curve defined by -x^3 - x*z^2 + y^2*z
```
"""
mutable struct ProjEllipticCurve{S} <: ProjectivePlaneCurve{S}
  eq::Oscar.MPolyElem_dec{S}
  degree::Int
  components::Dict{ProjEllipticCurve{S}, Int}
  point::Oscar.Geometry.ProjSpcElem{S}
  maps
  Hecke_ec::EllCrv{S}

  function ProjEllipticCurve{S}(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    ishomogeneous(eq) || error("The defining equation is not homogeneous")
    _iselliptic(eq) || error("Not an elliptic curve")
    isweierstrass_form(eq.f) || error("Not in Weierstrass form, please specify the point at infinity")
    v = shortformtest(eq.f)
    T = parent(eq)
    K = T.R.base_ring
    PP = proj_space(K, 2)
    V = gens(T)
    new{S}(eq, 3, Dict{ProjEllipticCurve{S}, Int}(), Oscar.Geometry.ProjSpcElem(PP[1], [K(0), K(1), K(0)]), [hom(T, T, V), hom(T, T, V)], Hecke.EllipticCurve(v[2], v[1]))
  end

  function ProjEllipticCurve(eq::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where {S <: FieldElem}
     nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
     iszero(evaluate(eq, P.v)) || error("The point is not on the curve")
     !isconstant(eq) || error("The defining equation must be non constant")
     ishomogeneous(eq) || error("The defining equation is not homogeneous")
     isinflection(eq, P) || error("Not an inflection point -- structure implemented only with an inflection point as base point.")
     _iselliptic(eq) || error("Not an elliptic curve")
     T = parent(eq)
     L = _change_coord(eq, P)
     H = L[1](eq)
     isweierstrass_form(H) || error("Not in Weierstrass form")
     v = shortformtest(H)
     new{S}(eq, 3, Dict{ProjEllipticCurve{S}, Int}(), P, L, Hecke.EllipticCurve(v[2], v[1]))
  end

  function ProjEllipticCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: Nemo.fmpz_mod}
    nvars(parent(eq)) == 3 || error("The defining equation must belong to a ring with three variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    ishomogeneous(eq) || error("The defining equation is not homogeneous")
    isweierstrass_form(eq.f) || error("Not in Weierstrass form")
    v = shortformtest(eq.f)
    v[1] || error("Not in short Weierstrass form")
    d = _discr(eq)
    T = parent(eq)
    K = T.R.base_ring
    n = modulus(K)
    gcd(Hecke.data(d), n) == ZZ(1) || error("The discriminant is not invertible")
    PP = proj_space(K, 2)
    V = gens(T)
    E = new{S}()
    E.eq = eq
    E.degree = 3
    E.point = Oscar.Geometry.ProjSpcElem(PP[1], [K(0), K(1), K(0)])
    E.Hecke_ec = Hecke.EllipticCurve(v[2], v[1])
    return E
  end
end

ProjEllipticCurve(eq::Oscar.MPolyElem_dec{S}) where {S <: FieldElem} = ProjEllipticCurve{S}(eq)

function ProjEllipticCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
  R, _ = grade(parent(eq))
  return ProjEllipticCurve{S}(R(eq))
end

function Base.show(io::IO, C::ProjEllipticCurve)
  if !get(io, :compact, false)
     print(io, "Projective elliptic curve defined by ", C.eq)
  else
     print(io, C.eq)
  end
end

################################################################################
################################################################################

function curve_components(E::ProjEllipticCurve{S}) where S <: FieldElem
   E.components = Dict(E => 1)
   return E.components
end

################################################################################

@doc Markdown.doc"""
    weierstrass_form(E::ProjEllipticCurve{S}) where {S <: FieldElem}

Return the equation of a projective elliptic curve defined by an equation in
Weierstrass form and which is linearly equivalent to `E`.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> F = T(-x^3 - 3*x^2*y - 3*x*y^2 - x*z^2 - y^3 + y^2*z - y*z^2 - 4*z^3)
-x^3 - 3*x^2*y - 3*x*y^2 - x*z^2 - y^3 + y^2*z - y*z^2 - 4*z^3

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
(-1 : 1 : 0)

julia> E = Oscar.ProjEllipticCurve(F, P)
Projective elliptic curve defined by -x^3 - 3*x^2*y - 3*x*y^2 - x*z^2 - y^3 + y^2*z - y*z^2 - 4*z^3

julia> Oscar.weierstrass_form(E)
-x^3 - x*z^2 + y^2*z - 4*z^3
```
"""
function weierstrass_form(E::ProjEllipticCurve{S}) where S <: FieldElem
   V = E.Hecke_ec.coeff
   R = parent(E.eq)
   x = gen(R, 1)
   y = gen(R, 2)
   z = gen(R, 3)
   if length(V) == 5
      F = R(y^2*z + V[1]*x*y*z + V[3]*y*z^2 - x^3 - V[2]*x^2*z - V[4]*x*z^2 - V[5]*z^3)
   elseif length(V) == 2
      F = R(y^2*z - x^3 - V[1]*x*z^2 - V[2]*z^3)
   end
   return F
end


function base_point(E::ProjEllipticCurve)
   return E.point
end

################################################################################
################################################################################

@doc Markdown.doc"""
    discriminant(E::ProjEllipticCurve{S}) where S <: FieldElem

Return the discriminant of the projective elliptic curve `E`.
"""
function Oscar.discriminant(E::ProjEllipticCurve{S}) where S <: FieldElem
   return Hecke.discriminant(E.Hecke_ec)
end

function Oscar.issmooth(E::ProjEllipticCurve{S}) where {S <: FieldElem}
   return true
end

@doc Markdown.doc"""
    j_invariant(E::ProjEllipticCurve{S}) where S <: FieldElem

Return the j-invariant of the projective elliptic curve `E`.
"""
function Oscar.j_invariant(E::ProjEllipticCurve{S}) where S <: FieldElem
   return j_invariant(E.Hecke_ec)
end

################################################################################
# Structure for points on elliptic curve.
################################################################################
#Point_EllCurve(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where {S <: Nemo.fmpz_mod}

@doc Markdown.doc"""
    Point_EllCurve(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where {S <: FieldElem}
    Point_EllCurve(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where {S <: Nemo.fmpz_mod}

Create the point `P` on the elliptic curve `E`.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> E = Oscar.ProjEllipticCurve(T(y^2*z - x^3 + 2*x*z^2))
Projective elliptic curve defined by -x^3 + 2*x*z^2 + y^2*z

julia> PP = Oscar.proj_space(E)
Projective space of dim 2 over Rational Field

julia> P = Oscar.Geometry.ProjSpcElem(PP, [QQ(2), QQ(2), QQ(1)])
(2 : 2 : 1)

julia> Q = Oscar.Point_EllCurve(E, P)
(2 : 2 : 1)
```
"""
struct Point_EllCurve{S <: RingElem}
   Pt::Oscar.Geometry.ProjSpcElem{S}
   C::ProjEllipticCurve{S}
   Hecke_Pt::Hecke.EllCrvPt{S}

   function Point_EllCurve(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S
      P in E || error("The point is not on the curve")
      new{S}(P, E, _point_toweierstrass(E, P))
   end

   function Point_EllCurve(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where {S <: Nemo.fmpz_mod}
      T = E.eq.parent
      A = base_ring(T)
      zz = gen(T, 3)
      P in E || error("The point is not on the curve")
      if iszero(P.v[1]) && iszero(P.v[3])
         Q = new{S}(P, E, Hecke.infinity(E.Hecke_ec))
      else
         P.v[3] == A(1) || error("Requires a point satisfying ", zz ," = 1")
         Q = new{S}(P, E, _point_tohecke_znz(E, P))
      end
      return Q
   end
end

function Base.show(io::IO, P::Point_EllCurve)
   print(io, P.Pt)
end

################################################################################

function _image_point(phi, V::Vector{S}) where S <: FieldElem
   return [evaluate(Oscar._images(phi)[i], V) for i in 1:3]
end

################################################################################

function _point_toweierstrass(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
   P in E || error("not on the curve")
   L = E.maps
   V = _image_point(L[2], P.v)
   if iszero(V[3])
      return Hecke.infinity(E.Hecke_ec)
   else
      W = [V[1]//V[3], V[2]//V[3]]
      return E.Hecke_ec(W)
   end
end

function _point_tohecke_znz(E::ProjEllipticCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: Nemo.fmpz_mod
   if iszero(P.v[3])
      return Hecke.infinity(E.Hecke_ec)
   else
      return E.Hecke_ec([P.v[1], P.v[2]])
   end
end

################################################################################

function _point_fromweierstrass(E::ProjEllipticCurve{S}, PP::Oscar.Geometry.ProjSpc{S}, P::Hecke.EllCrvPt{S}) where S <: FieldElem
   E.Hecke_ec.coeff == P.parent.coeff || error("not the same curve")
   L = E.maps
   K = P.parent.base_field
   if P.isinfinite
      V = [K(0), K(1), K(0)]
   else
      V = [P.coordx, P.coordy, K(1)]
   end
   W = _image_point(L[1], V)
   return Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, W))
end

################################################################################

function _EllCrvPt(P::Point_EllCurve{S}) where S <: FieldElem
   return P.Hecke_Pt
end

################################################################################

@doc Markdown.doc"""
    curve(P::Point_EllCurve{S}) where S <: FieldElem

Return the curve on which the point `P` is considered.
"""
function curve(P::Point_EllCurve{S}) where S <: FieldElem
   return P.C
end

################################################################################
@doc Markdown.doc"""
    proj_space(P::Point_EllCurve{S}) where S <: FieldElem

Return the projective space to which the point `P` belongs.
"""
function Oscar.Geometry.proj_space(P::Point_EllCurve{S}) where S <: FieldElem
   return P.Pt.parent
end

################################################################################

function ==(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: FieldElem
   return P.C == Q.C && P.Pt == Q.Pt
end

################################################################################
function Base.:+(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: FieldElem
   E = P.C
   PP = proj_space(P)
   A = _EllCrvPt(P) + _EllCrvPt(Q)
   return _point_fromweierstrass(E, PP, A)
end

################################################################################

function Base.:-(P::Point_EllCurve{S}) where S <: FieldElem
   E = P.C
   PP = proj_space(P)
   return _point_fromweierstrass(E, PP, -_EllCrvPt(P))
end

function Base.:-(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: FieldElem
   E = P.C
   PP = proj_space(P)
   A = _EllCrvPt(P) + (-_EllCrvPt(Q))
   return _point_fromweierstrass(E, PP, A)
end
################################################################################

function Base.:*(n::Int64, P::Point_EllCurve{S}) where S <: FieldElem
   E = curve(P)
   PP = proj_space(P)
   return _point_fromweierstrass(E, PP, n*_EllCrvPt(P))
end

################################################################################
# helping function

function contract(i::Int, F::Oscar.MPolyElem{S}) where {S <: FieldElem}
   R = parent(F)
   P = MPolyBuildCtx(R)
   L = exponent_vectors(F)
   e = zeros(Int64, nvars(R)-1)
   insert!(e, i, -1)
   A = [v + e for v in L]
   Co = coefficients(F)
   B = [c for c in Co]
   C = []
   E = []
   j = 0
   for v in A
      j = j + 1
      if v[i] >= 0
         push!(C, B[j])
         push!(E, v)
      end
   end
   Z = zip(C, E)
   for (c, v) in Z
      push_term!(P, c, v)
   end
   return finish(P)
end

################################################################################
# helping function

function _eva(F::Oscar.MPolyElem{S}) where {S <: FieldElem}
   R = parent(F)
   nvars(R) == 3 || error("wrong number of variables")
   X = gens(R)
   return evaluate(evaluate(F, [X[1], X[2], R(0)]), [R(1), X[1], R(1)])
end

################################################################################
# This code is inspired from Macaulay2, based on Tibouchi, M. ''A Nagell Algorithm
# in Any Characteristic'', Cryptography and Security: From Theory to
# Applications, 474--479, 2012.

@doc Markdown.doc"""
    toweierstrass(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Given a smooth plane cubic projective curve `C` and a point `P` on the curve,
return an elliptic curve birationally equivalent to `C` given by an equation in
long Weierstrass form.

# Examples
```jldoctest
julia> S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> PP = proj_space(QQ, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

julia> Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
(-1 : 1 : 0)

julia> D = Oscar.ProjPlaneCurve(T(-x^3 - 3*x^2*y + 2*x^2*z - 3*x*y^2 + 3*x*y*z - 4*x*z^2 - y^3 - y*z^2 + 6*z^3))
Projective plane curve defined by -x^3 - 3*x^2*y + 2*x^2*z - 3*x*y^2 + 3*x*y*z - 4*x*z^2 - y^3 - y*z^2 + 6*z^3

julia> Oscar.toweierstrass(D, Q)
-x^3 - 2*x^2*z + x*y*z - 4*x*z^2 + y^2*z + 3*y*z^2 - 6*z^3
```
"""
function toweierstrass(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
   _iselliptic(C.eq) || error("not elliptic curve")
   P in C || error("The point ", P, " is not on the curve")
   F = defining_equation(C)
   R = parent(F)
   T = parent(C.eq)
   X = gens(R)
   PP = P.parent
   if iszero(P.v[2])
      if iszero(P.v[1])
         F = evaluate(F, [X[1], X[3], X[2]])
         Q = Oscar.Geometry.ProjSpcElem(PP, [P.v[1], P.v[3], P.v[2]])
      else
         F = evaluate(F, [X[2], X[1], X[3]])
         Q = Oscar.Geometry.ProjSpcElem(PP, [P.v[2], P.v[1], P.v[3]])
      end
   else
      Q = P
   end
   F = evaluate(F, [Q.v[2]*X[1] + Q.v[1]*X[2], Q.v[2]*X[2], Q.v[2]*X[3] + Q.v[3]*X[2]])
   # coeff of x*y^2, y^2*z
   a = coeff(F, X[1]*X[2]^2)
   b = coeff(F, X[2]^2*X[3])
   if iszero(a)
      F = evaluate(F, [X[3], X[2], X[1]])
   else
      F = evaluate(F, [a*X[1] - b*X[3], a*X[2], a*X[3]])
   end
   # coeff of y*z^2, z^3
   c = coeff(F, X[2]*X[3]^2)
   d = coeff(F, X[3]^3)
   if iszero(c)
      F = evaluate(F, [R(1), X[2], X[1]])
      F = F*1//coeff(F, X[2]^2)
      G = [evaluate(contract(2, F), [X[1], R(0), X[3]])]
      G = push!(G, -evaluate(F, [X[1], R(0), X[3]]))
   else
      F = evaluate(F, [c*X[1], c*X[2] - d*X[3], c*X[3]])
      ff = [_eva(contract(3, contract(3, F)))]
      ff = push!(ff, _eva(contract(3, F)))
      ff = push!(ff, _eva(F))
      G = [ff[2]]
      G = push!(G, -ff[1]*ff[3])
   end
   u = leading_coefficient(G[2])
   a1 = evaluate(contract(1, G[1]), [R(0), X[2], X[3]])
   a3 = u*evaluate(G[1], [R(0), X[2], X[3]])
   a2 = evaluate(contract(1, contract(1, G[2])), [R(0), X[2], X[3]])
   a4 = u*evaluate(contract(1, G[2]), [R(0), X[2], X[3]])
   a6 = u^2*evaluate(G[2], [R(0), X[2], X[3]])
   return T(X[2]^2*X[3] + a1*X[1]*X[2]*X[3] + a3*X[2]*X[3]^2 - X[1]^3 - a2*X[1]^2*X[3] - a4*X[1]*X[3]^2 - a6*X[3]^3)
end


################################################################################
@doc Markdown.doc"""
    proj_space(E::ProjEllipticCurve{S}) where S <: FieldElem

Return the projective space to which the base point of the elliptic curve `E` belongs.
"""
function Oscar.Geometry.proj_space(E::ProjEllipticCurve{S}) where S <: FieldElem
   return base_point(E).parent
end

################################################################################
#
# Interface calling the functions in Hecke
#
################################################################################

@doc Markdown.doc"""
    order(E::ProjEllipticCurve{S}) where S <: FieldElem

Given an elliptic curve `E` over a finite field $\mathbf F$, computes
$\#E(\mathbf F)$.
"""
function Oscar.order(E::ProjEllipticCurve{S}) where S <: FieldElem
   return Hecke.order(E.Hecke_ec)
end

################################################################################

@doc Markdown.doc"""
    rand(E::ProjEllipticCurve{S}) where S <: FieldElem

Return a random point on the elliptic curve `E` defined over a finite field.
"""
function Oscar.rand(E::ProjEllipticCurve{S}) where S <: FieldElem
   PP = proj_space(E)
   H = E.Hecke_ec
   if !Hecke.isshort(H)
      L = Hecke.short_weierstrass_model(H)
      P = Hecke.rand(L[1])
      Q = L[3](P)
   else
      Q = Hecke.rand(H)
   end
   return _point_fromweierstrass(E, PP, Q)
end

################################################################################

@doc Markdown.doc"""
    list_rand(E::ProjEllipticCurve, N::Int)

Return a list of `N` random points on the elliptic curve `E` defined over a finite field.
"""
function list_rand(E::ProjEllipticCurve, N::Int)
   PP = proj_space(E)
   H = E.Hecke_ec
   M = []
   if !Hecke.isshort(H)
      L = Hecke.short_weierstrass_model(H)
      for i in 1:N
         P = Hecke.rand(L[1])
         Q = L[3](P)
         push!(M, Q)
      end
   else
      for i in 1:N
         Q = Hecke.rand(H)
         push!(M, Q)
      end
   end
   return [_point_fromweierstrass(E, PP, Q) for Q in M]
end

################################################################################
@doc Markdown.doc"""
    order(P::Point_EllCurve{fmpq})

Return the order of the point `P` or `0` if the order is infinite.
"""
function Oscar.order(P::Point_EllCurve{fmpq})
   return Hecke.order(P.Hecke_Pt)
end

################################################################################

@doc Markdown.doc"""
    istorsion_point(P::Point_EllCurve{fmpq})

Return whether the point `P` is a torsion point.
"""
function Oscar.istorsion_point(P::Point_EllCurve{fmpq})
   return Hecke.istorsion_point(P.Hecke_Pt)
end

################################################################################
@doc Markdown.doc"""
    torsion_points_lutz_nagell(E::ProjEllipticCurve{fmpq})

Return the rational torsion points of the elliptic curve `E` using the
Lutz-Nagell theorem.
"""
function Oscar.torsion_points_lutz_nagell(E::ProjEllipticCurve{fmpq})
   PP = proj_space(E)
   H = E.Hecke_ec
   M = Hecke.torsion_points_lutz_nagell(H)
   return [_point_fromweierstrass(E, PP, a) for a in M]
end

################################################################################

@doc Markdown.doc"""
    torsion_points_division_poly(E::ProjEllipticCurve{fmpq})

Return the rational torsion points of a rational elliptic curve `E` using
division polynomials.
"""
function Oscar.torsion_points_division_poly(E::ProjEllipticCurve{fmpq})
   PP = proj_space(E)
   H = E.Hecke_ec
   M = Hecke.torsion_points_division_poly(H)
   return [_point_fromweierstrass(E, PP, a) for a in M]
end
