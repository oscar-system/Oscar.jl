```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["divisors.md"]
```

# Elliptic Curves

An elliptic plane curve is a projective plane curve of degree 3 together with a point of the curve, called the base point. An operation of addition of points can be defined on elliptic curves.

## Constructor

An elliptic curve is a subtype of the abstract type ProjectivePlaneCurve. To define an elliptic curve over a field, one can either give as input an equation and the point at infinity, or just an equation in Weierstrass form. In the latter case, the point at infinity is ``(0 : 1 : 0)``.

Considering elliptic curves over a ring is helpful in some primality proving test. We introduce here a structure of elliptic curve over a ring. In that case, we always assume the equation to be in Weierstrass form, with infinity point ``(0 : 1 : 0)``.

```@docs
ProjEllipticCurve
```

## Points on Elliptic Curves

We define a specific structure for the points on an elliptic curve, on which the operation of addition and multiplication by an integer are defined.

```@docs
Point_EllCurve
```

## Methods

Most of the functions described for projective plane curves are also available for elliptic curves over a field. We describe here the functions which are specific to elliptic curves.

### Basic functions

```@docs
proj_space(E::ProjEllipticCurve{S}) where S <: FieldElem
proj_space(P::Point_EllCurve{S}) where S <: FieldElem
curve(P::Point_EllCurve{S}) where S <: FieldElem
```

Addition and multiplication by an integer of a point on an elliptic curve can be performed using the usual symbols `+` and `*`.

#### Example

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

julia> Q+Q
(9//4 : -21//8 : 1)

julia> 3*Q
(338 : 6214 : 1)

```

### Weierstrass form

```@docs
weierstrass_form(E::ProjEllipticCurve{S}) where {S <: FieldElem}
toweierstrass(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

### Invariant of Elliptic Curves, torsion points...

```@docs
discriminant(E::ProjEllipticCurve{S}) where S <: FieldElem
j_invariant(E::ProjEllipticCurve{S}) where S <: FieldElem
is_torsion_point(P::Point_EllCurve{fmpq})
torsion_points_lutz_nagell(E::ProjEllipticCurve{fmpq})
torsion_points_division_poly(E::ProjEllipticCurve{fmpq})
order(P::Point_EllCurve{fmpq})
```

The following functions are implemented for elliptic curves over finite fields:

```@docs
rand(E::ProjEllipticCurve{S}) where S <: FieldElem
list_rand(E::ProjEllipticCurve, N::Int)
order(E::ProjEllipticCurve{S}) where S <: FieldElem
```

The addition of points is not well defined for projective elliptic curves over a ring, that's why this case has to be considered separately. The following methods can for example be used for teaching purposes to show the steps of the Elliptic Curve Method.

```@docs
sum_Point_EllCurveZnZ(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: Nemo.fmpz_mod
IntMult_Point_EllCurveZnZ(m::fmpz, P::Point_EllCurve{S}) where S <: Nemo.fmpz_mod
rand_pair_EllCurve_Point(R::Oscar.MPolyRing_dec{S}, PP::Oscar.Geometry.ProjSpc{S}) where S <: Nemo.fmpz_mod
```

### Primality proving

#### Elliptic Curve Method

Projective elliptic curves over a ring are for example used in the Elliptic Curve Method. We give here an example (see Example 7.1 of [Was08](@cite)) on how to use the previous functions to apply it. 

```jldoctest
julia> n = 4453
4453

julia> A = ResidueRing(ZZ, ZZ(n))
Integers modulo 4453

julia> S, (x,y,z) = PolynomialRing(A, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integers modulo 4453, AbstractAlgebra.Generic.MPoly{fmpz_mod}[x, y, z])

julia> T, _ = grade(S)
(Multivariate Polynomial Ring in x, y, z over Integers modulo 4453 graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpz_mod, AbstractAlgebra.Generic.MPoly{fmpz_mod}}[x, y, z])

julia> E = Oscar.ProjEllipticCurve(T(y^2*z - x^3 - 10*x*z^2 + 2*z^3))
Projective elliptic curve defined by 4452*x^3 + 4443*x*z^2 + y^2*z + 2*z^3

julia> PP = proj_space(A, 2)
(Projective space of dim 2 over Integers modulo 4453
, MPolyElem_dec{fmpz_mod, AbstractAlgebra.Generic.MPoly{fmpz_mod}}[x[0], x[1], x[2]])

julia> Q = Oscar.Geometry.ProjSpcElem(PP[1], [A(1), A(3), A(1)])
(1 : 3 : 1)

julia> P = Oscar.Point_EllCurve(E, Q)
(1 : 3 : 1)

julia> P2 = Oscar.IntMult_Point_EllCurveZnZ(ZZ(2), P)
(4332 : 3230 : 1)

julia> Oscar.sum_Point_EllCurveZnZ(P, P2)
ERROR: 61 is not invertible in the base ring, cannot perform the sum

```

The last sum is not defined, and the error which is shown when we ask for the sum gives us a factor of `4453`. The Elliptic Curve Method is implemented and can be called using:

```@docs
ECM(n::fmpz; nbcurve::Int = 25000, multfact::fmpz = factorial(ZZ(10^4)))
```

#### Elliptic Curve Primality Proving

Elliptic curves over finite fields or rings are used for the method called "Elliptic Curve Primality Proving". We implemented here the version relying on Atkin-Morain's test. The present implementation is not intended to be competitive.

```@docs
ECPP(n::fmpz)
```

#### Other algorithms

TODO: The following algorithms are not directly related to Plane Curves, they might be moved to another section of the documentation.

```@docs
cornacchia_algorithm(d::fmpz, m::fmpz)
Miller_Rabin_test(N::fmpz, k::Int64 = 20)
Pollard_rho(N::fmpz, bound::Int = 50000)
Pollard_p_1(N::fmpz, B::fmpz = ZZ(10)^5)
```



