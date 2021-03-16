export AffineEllipticCurve, discriminant, issmooth, j_invariant

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
   return Hecke.disc(E.Hecke_ec)
end

function Oscar.issmooth(E::AffineEllipticCurve)
   return true
end

function Oscar.j_invariant(E::AffineEllipticCurve)
   return j_invariant(E.Hecke_ec)
end
