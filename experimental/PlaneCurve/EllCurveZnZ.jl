export sum_Point_EllCurveZnZ, ECM, rand_pair_EllCurve_Point,
       IntMult_Point_EllCurveZnZ

################################################################################
# Elliptic curves over a ring Z/nZ
################################################################################
################################################################################
# Helping functions
################################################################################
# addition on two projective points on an elliptic curve over Z/nZ, given in
# short Weierstrass form, when it is possible. It gives as an output a pair
# composed of the coordinates of the sum and ZZ(1) if the sum exists, and
# (0, 0, 0) and the gcd otherwise.
# Adapted from the corresponding code in Hecke.

function _add(P::Array{Nemo.fmpz_mod, 1}, Q::Array{Nemo.fmpz_mod, 1}, E::Array{Nemo.fmpz_mod, 1})
   length(P) == 3 && length(Q) == 3 || error("arrays of size 3 required")
   length(E) == 2 || error("array of size 2 required")
   A = P[1].parent
   A == Q[1].parent && A == E[1].parent || error("Not the same parent")
   P[3] != A(1) && P != [A(0), A(1), A(0)] && error("require infinity point or last coordinate 1")
   Q[3] != A(1) && Q != [A(0), A(1), A(0)] && error("require infinity point or last coordinate 1")
   n = modulus(A)
   if P == [A(0), A(1), A(0)]
      return [Q, ZZ(1)]
   elseif Q == [A(0), A(1), A(0)]
      return [P, ZZ(1)]
   end
   if P[1] != Q[1]
      d = Q[1] - P[1]
      g = gcd(n, Hecke.data(d))
      if g == ZZ(1)
         m = divexact(Q[2] - P[2], d)
         x = m^2 - P[1] - Q[1]
         y = m*(P[1] - x) - P[2]
         z = A(1)
      else
         x = A(0)
         y = A(0)
         z = A(0)
      end
   elseif P[2] != Q[2]
      return [[A(0), A(1), A(0)], ZZ(1)]
   elseif P[2] != 0
      d = 2*P[2]
      g = gcd(n, Hecke.data(d))
      if g == ZZ(1)
         m = divexact(3*(P[1])^2 + E[1], 2*P[2])
         x = m^2 - 2*P[1]
         y = m*(P[1] - x) - P[2]
         z = A(1)
      else
         x = A(0)
         y = A(0)
         z = A(0)
      end
   else
      return [[A(0), A(1), A(0)], ZZ(1)]
   end
   return [[x, y, z], g]
end

################################################################################
# Creates a pair composed of the coordinates of a point and the coefficients of
# an elliptic curve going through that point.

function _rand_point_curve(A::Nemo.FmpzModRing)
   n = modulus(A)
   a = A(rand(ZZ(1):n))
   P = [A(rand(ZZ(1):n)), A(rand(ZZ(1):n)), A(ZZ(1))]
   b = A(P[2]^2 - P[1]^3 - a*P[1])
   return [P, [a, b]]
end

################################################################################
# Creates a list of pairs composed of a the coefficitents of an elliptic curve
# and the coordinates of a point on it.

function _rand_list(A::Nemo.FmpzModRing, N::Int64)
   E = Array{Nemo.fmpz_mod, 1}[]
   P = Array{Nemo.fmpz_mod, 1}[]
   for i in 1:N
      L = _rand_point_curve(A)
      E = push!(E, L[2])
      P = push!(P, L[1])
   end
   return [E, P]
end

################################################################################
# Returns the coordinates of m*P and ZZ(1) when the computation is possible, and
# returns (0, 0, 0) and the gcd otherwise.

function _scalar_mult(P::Array{Nemo.fmpz_mod, 1}, E::Array{Nemo.fmpz_mod, 1}, m::Int64)
   length(P) == 3 || error("arrays of size 3 required")
   length(E) == 2 || error("array of size 2 required")
   A = P[1].parent
   if m < 0
      error("Positive integer expected")
   elseif m == 0
      return [[A(0), A(1), A(0)], ZZ(1)]
   elseif m == 1
      return [P, ZZ(1)]
   else
      Q = _add(P, P, E)
      if Q[2] != ZZ(1)
         return Q
      elseif iszero(m-2)
         return Q
      else
         for i in 1:m-1
            Q = _add(Q[1], P, E)
            if Q[2] != ZZ(1)
               return Q
            end
         end
      end
      return Q
   end
end

################################################################################
# Returns the coordinates of m!*P and ZZ(1) when the computation is possible, and
# returns (0, 0, 0) and the gcd otherwise.

function _fac_mult(P::Array{Nemo.fmpz_mod, 1}, E::Array{Nemo.fmpz_mod, 1}, m::Int64)
   length(P) == 3 || error("arrays of size 3 required")
   length(E) == 2 || error("array of size 2 required")
   A = P[1].parent
   if m < 0
      error("Positive integer expected")
   elseif m == 1 || m == 0
      return [P, ZZ(1)]
   else
      Q = _scalar_mult(P, E, 2)
      if Q[2] != ZZ(1)
         return Q
      elseif iszero(m-2)
         return Q
      else
         for i in 1:m-2
            Q = _scalar_mult(Q[1], E, i+2)
            if Q[2] != ZZ(1)
               return Q
            end
         end
      end
   end
   return Q
end

################################################################################
# Returns the elliptic plane curve corresponding to the array.

function _toProjEllipticCurve(R::MPolyRing{S}, E::Array{S, 1}) where S <: Nemo.fmpz_mod
   length(E) == 2 || error("array of size 2 required")
   x = gen(R, 1)
   y = gen(R, 2)
   z = gen(R, 3)
   return ProjEllipticCurve(R(y^2*z - x^3 - E[1]*x*z^2 - E[2]*z^3))
end

################################################################################
# Functions
################################################################################
@doc Markdown.doc"""
    ECM(n::Int64; nbcurve::Int64 = 10, multfact::Int64 = 8)

Return a factor of `n`, obtained with the Elliptic Curve Method.
"""
function ECM(n::fmpz; nbcurve::Int64 = 10, multfact::Int64 = 8)
   A = ResidueRing(ZZ, n)
   N = nbcurve
   m = multfact
   L = _rand_list(A, N)
   for i in 1:length(L)
      Q = _fac_mult(L[2][i], L[1][i], m)
      if Q[2] != ZZ(1)
         return Q[2]
      end
   end
   return ZZ(1)
end

################################################################################
# the sum of two points might not be a point, this is not a group operation.

@doc Markdown.doc"""
    sum_Point_EllCurveZnZ(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: Nemo.fmpz_mod

Return, if possible, the sum of the points `P` and `Q`, and an error otherwise.
"""
function sum_Point_EllCurveZnZ(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: Nemo.fmpz_mod
   A = P.Pt[1].parent
   E = P.C
   E.Hecke_ec.short || error("requires short Weierstrass form")
   Q = _add(P.Pt.v, Q.Pt.v, E.Hecke_ec.coeff)
   PP = P.Pt.parent
   if Q[2] == A(1)
      return Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, Q[1]))
   else
      error(Q[2], " is not invertible in the base ring, cannot perform the sum")
   end
end

################################################################################
@doc Markdown.doc"""
    IntMult_Point_EllCurveZnZ(m::Int, P::Point_EllCurve{S}) where S <: Nemo.fmpz_mod

Return, if possible, the point `mP`, and an error otherwise.
"""
function IntMult_Point_EllCurveZnZ(m::Int, P::Point_EllCurve{S}) where S <: Nemo.fmpz_mod
   E = P.C
   E.Hecke_ec.short || error("requires short Weierstrass form")
   PP = P.Pt.parent
   Q = _scalar_mult(P.Pt.v, E.Hecke_ec.coeff, m)
   if Q[2] == ZZ(1)
      return Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, Q[1]))
   else
      error(Q[2], " is not invertible in the base ring, cannot perform the sum")
   end
end

################################################################################
@doc Markdown.doc"""
    rand_pair_ellcurve_point(R::Oscar.MPolyRing_dec{S}, PP::Oscar.Geometry.ProjSpc{S}) where S <: Nemo.fmpz_mod

Return a pair composed of an elliptic plane curve `E` with equation in `R`,
and a point `P` on `E`.
"""
function rand_pair_EllCurve_Point(R::Oscar.MPolyRing_dec{S}, PP::Oscar.Geometry.ProjSpc{S}) where S <: Nemo.fmpz_mod
   A = base_ring(R)
   n = modulus(A)
   i = 0
   L = _rand_point_curve(A)
   while i < 100 && !isunit(4*L[2][1]^3 + 27*L[2][2]^2)
      i = i+1
      L = _rand_list(A, 1)
   end
   if isunit(4*L[2][1]^3 + 27*L[2][2]^2)
      E = _toProjEllipticCurve(R, L[2])
      return [E, Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, L[1]))]
   else
      error("Did not manage do find an elliptic curve with invertible discriminant")
   end
end
