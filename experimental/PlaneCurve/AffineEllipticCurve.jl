export AffineEllipticCurve, discriminant, issmooth, j_invariant,
       iselliptic, toweierstrass

################################################################################
# Helping functions
################################################################################

function isweierstrassform(eq::Oscar.MPolyElem{S}) where S <: FieldElem
   R = parent(eq)
   x = gen(R, 1)
   y = gen(R, 2)
   A = [x*y, y, x^2, x, R(1)]
   f = deepcopy(eq)
   f = f - y^2 + x^3
   for t in A
      f = f - coeff(eq, t)*t
   end
   return iszero(f)
end

function shortformtest(eq::Oscar.MPolyElem{S}) where S <: FieldElem
   R = parent(eq)
   x = gen(R, 1)
   y = gen(R, 2)
   A = [coeff(eq, x*y), -coeff(eq, x^2), coeff(eq, y), -coeff(eq, x), -coeff(eq, R(1))]
   iszero(A[1]) && iszero(A[2]) && iszero(A[3]) && return (true, [A[4], A[5]])
   return (false, A)
end

################################################################################
# Definition EllipticCurve
################################################################################

mutable struct AffineEllipticCurve{S <: FieldElem} <: PlaneCurve
  eq::Oscar.MPolyElem{S}                # Equation of the curve (polynomial in two variables)
  degree::Int                           # degree of the equation of the curve
  components::Dict{AffineEllipticCurve{S}, Int}
  Hecke_ec::EllCrv{S}
  function AffineEllipticCurve{S}(eq::Oscar.MPolyElem{S}) where {S <: FieldElem}
    nvars(parent(eq)) == 2 || error("The defining equation must belong to a ring with two variables")
    !isconstant(eq) || error("The defining equation must be non constant")
    !isweierstrassform(eq) && error("Not in Weierstrass form!")
    v = shortformtest(eq)
    new{S}(eq,
           -1,                   # -1 when it is not computed yet
            Dict{AffineEllipticCurve{S}, Int}(),
            Hecke.EllipticCurve(v[2], v[1]))
  end
end

AffineEllipticCurve(eq::Oscar.MPolyElem{S}) where {S <: FieldElem} = AffineEllipticCurve{S}(eq)

function Base.show(io::IO, C::AffineEllipticCurve)
  if !get(io, :compact, false)
     println(io, "Affine elliptic curve defined by ", C.eq)
  else
     println(io, C.eq)
  end
end

################################################################################
################################################################################

function Oscar.discriminant(E::AffineEllipticCurve)
   return Hecke.discriminant(E.Hecke_ec)
end

function Oscar.issmooth(E::AffineEllipticCurve)
   return true
end

function Oscar.j_invariant(E::AffineEllipticCurve)
   return j_invariant(E.Hecke_ec)
end

################################################################################

function iselliptic(C::PlaneCurve)
   return degree(C) == 3 && issmooth_curve(C)
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
   Co = coeffs(F)
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
# This code is copied from Macaulay2, based on Tibouchi, M. ''A Nagell Algorithm
# in Any Characteristic'', Cryptography and Security: From Theory to
# Applications, 474--479, 2012.

@doc Markdown.doc"""
    toweierstrass(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem

Given a smooth plane cubic projective curve `C` and a point `P` on the curve,
return an equation in long Weierstrass form of an elliptic curve isormorphic to
 `C`.
"""
function toweierstrass(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
   iselliptic(C) || error("not elliptic curve")
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
   u = lc(G[2])
   a1 = evaluate(contract(1, G[1]), [R(0), X[2], X[3]])
   a3 = u*evaluate(G[1], [R(0), X[2], X[3]])
   a2 = evaluate(contract(1, contract(1, G[2])), [R(0), X[2], X[3]])
   a4 = u*evaluate(contract(1, G[2]), [R(0), X[2], X[3]])
   a6 = u^2*evaluate(G[2], [R(0), X[2], X[3]])
   return T(X[2]^2*X[3] + a1*X[1]*X[2]*X[3] + a3*X[2]*X[3]^2 - X[1]^3 - a2*X[1]^2*X[3] - a4*X[1]*X[3]^2 - a6*X[3]^3)
end
