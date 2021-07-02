```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["exp_plane_curves.md"]
```

# Plane Curves

Affine plane curves and projective plane curves are considered in this module. An affine plane curve is defined by a polynomial in two variables, whereas a projective plane curve is defined by an homogenous polynomial belonging to a graded polynomial ring in three variables.  

## Points

Several functions involve a point on a curve. We describe in this section how to deal with points, both in the affine and projective settings.

### Points in the affine space

A point in the affine space can be defined as follows:
```@docs
Point
```
We consider also the following function for points.

```@docs
ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem
```

### Points in the projective space

In order to define a point in the projective plane, one needs first to define the projective plane as follows, where `K` is the base ring:

#### Example
```@repl oscar
K = QQ
PP = projective_space(K, 2)
```

Then, one can define a projective point as follows:

#### Example
```@repl oscar
P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(2), QQ(-5)])
```
## Affine and Projective Plane Curves

We consider two kinds of plane curves : affine and projective.

### Constructors

An affine plane curve is defined as the class of a two-variables polynomial ``C`` over a field ``K``, modulo the equivalence relation ``C \sim D \iff \exists \lambda \in K\backslash \{0\}, C = \lambda \cdot D``.  

#### Example

```@docs
AffinePlaneCurve
```

Similarly, a projective plane curve is defined as the class of a three-variables homogeneous polynomial ``C`` over a field ``K``, modulo the equivalence relation ``C\sim D \iff \exists \lambda \in K\backslash \{0\}, C = \lambda \cdot D``. The defining equation is supposed to belong to a graded ring.

```@docs
ProjPlaneCurve
```

A particular kind of projective curves is the case of elliptic curves, see the corresponding section for more information.


### Methods

The following methods are implemented for both affine and projective plane curves.

```@docs
in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem
```
```@docs
in(P::Oscar.Geometry.ProjSpcElem{S}, C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
```


```@docs
degree(C::Oscar.PlaneCurveModule.PlaneCurve)
```

```@docs
jacobi_ideal(C::Oscar.PlaneCurveModule.PlaneCurve)
```

```@docs
curve_components
```

```@docs
reduction
```

```@docs
isirreducible(C::Oscar.PlaneCurveModule.PlaneCurve{S}) where S <: FieldElem
```

```@docs
isreduced(C::Oscar.PlaneCurveModule.PlaneCurve{S}) where S <: FieldElem
```


```@docs
issmooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
```

```@docs
issmooth(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
```

```@docs
tangent(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
union(C::T, D::T) where T <: Oscar.PlaneCurveModule.PlaneCurve
```

```@docs
ring(C::Oscar.PlaneCurveModule.PlaneCurve)
```



```@docs
common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
```
```@docs
common_components(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
```

```@docs
curve_intersect(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
```

```@docs
curve_intersect([PP::Oscar.Geometry.ProjSpc{S}], C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
```

```@docs
intersection_multiplicity(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
```

```@docs
intersection_multiplicity(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
intersection_multiplicity
```

```@docs
curve_singular_locus(C::AffinePlaneCurve)
```

```@docs
curve_singular_locus([PP::Oscar.Geometry.ProjSpc{S}], C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
```

```@docs
multiplicity(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
```

```@docs
multiplicity(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
tangent_lines(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
```

```@docs
tangent_lines(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
aretransverse(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S<:FieldElem
```

```@docs
aretransverse(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S<:FieldElem
```


```@docs
issmooth_curve(C::AffinePlaneCurve)
```

```@docs
issmooth_curve(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve)
```

```@docs
arithmetic_genus(C::AffinePlaneCurve)
```

```@docs
arithmetic_genus(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve)
```

```@docs
geometric_genus(C::AffinePlaneCurve)
```

```@docs
geometric_genus(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
```

## Divisors

In order to consider divisors on curves, we restrict our attention to smooth and irreducible curve_components.

Let ``C`` be an affine or projective plane curve defined by an equation ``F``. Then any polynomial function ``G`` which is not divisible by ``F`` will vanish on ``C`` only at finitely many points. A way to encode these points together with their intersection multiplicities is to consider a divisor. A divisor on a curve is a formal finite sum of points of the curve with integer coefficients. A natural operation of addition can be defined on the set of divisors of a curve, which turns it into an Abelian group.


### Constructors

Divisors on curves are here introduced as a dictionary associating a point on the curve to its multiplicity.

```@docs
AffineCurveDivisor
```

```@docs
ProjCurveDivisor
```

### Methods

```@docs
curve(D::Oscar.PlaneCurveModule.CurveDivisor)
```

```@docs
degree(D::Oscar.PlaneCurveModule.CurveDivisor)
```

```@docs
iseffective(D::Oscar.PlaneCurveModule.CurveDivisor)
```

```@docs
curve_zero_divisor(C::ProjPlaneCurve{S}) where S <: FieldElem
```

```@docs
curve_zero_divisor(C::AffinePlaneCurve{S}) where S <: FieldElem
```

```@docs
multiplicity(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Point{S}) where {S <: FieldElem, T <: MPolyElem{S}}
```

```@docs
multiplicity(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}, P::Point{S}) where S <: FieldElem
```

```@docs
multiplicity(C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
multiplicity(C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Oscar.Geometry.ProjSpcElem{S})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
```

```@docs
divisor(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}) where S <: FieldElem
```

```@docs
divisor(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}) where {S <: FieldElem, T <: MPolyElem{S}}
```

```@docs
divisor([PP::Oscar.Geometry.ProjSpc{S}], C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}) where S <: FieldElem
```

```@docs
divisor(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
```


```@docs
islinearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
```

```@docs
isprincipal(D::ProjCurveDivisor{S}) where S <: FieldElem
```

```@docs
principal_divisor(D::ProjCurveDivisor{S}) where S <: FieldElem
```

```@docs
global_sections(D::ProjCurveDivisor)
```

```@docs
dimension_global_sections(D::ProjCurveDivisor)
```

## Elliptic Curves

An elliptic plane curve is a projective plane curve of degree 3 together with a point of the curve, called the base point. An operation of addition of points can be defined on elliptic curves.

### Constructors

An elliptic curve is a subtype of the abstract type ProjectivePlaneCurve. To define an elliptic curve over a field, one can either give as input an equation and the point at infinity, or just an equation in Weierstrass form. In the latter case, the point at infinity is ``(0 : 1 : 0)``.

Considering elliptic curves over a ring is helpful in some primality proving test. We introduce here a structure of elliptic curve over a ring. In that case, we always assume the equation to be in Weierstrass form, with infinity point ``(0 : 1 : 0)``.

```@docs
ProjEllipticCurve
```

We define a specific structure for the points on an elliptic curve, on which the operation of addition is definde.

```@docs
Point_EllCurve
```

### Methods

Most of the functions described for projective plane curves are also available for elliptic curves over a field. We describe here the functions which are specific to elliptic curves.

```@docs
weierstrass_form(E::ProjEllipticCurve{S}) where {S <: FieldElem}
```

```@docs
toweierstrass(C::ProjPlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
```

```@docs
discriminant(E::ProjEllipticCurve{S}) where S <: FieldElem
```

```@docs
j_invariant(E::ProjEllipticCurve{S}) where S <: FieldElem
```

```@docs
rand(E::ProjEllipticCurve{S}) where S <: FieldElem
```

```@docs
curve(P::Point_EllCurve{S}) where S <: FieldElem
```

```@docs
order(E::ProjEllipticCurve{S}) where S <: FieldElem
```


```@docs
istorsion_point(P::Point_EllCurve{fmpq})
```

```@docs
torsion_points_lutz_nagell(E::ProjEllipticCurve{fmpq})
```

```@docs
torsion_points_division_poly(E::ProjEllipticCurve{fmpq})
```

```@docs
order(P::Point_EllCurve{fmpq})
```

```@docs
list_rand(E::ProjEllipticCurve, N::Int)
```

```@docs
sum_Point_EllCurveZnZ(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: Nemo.fmpz_mod
```

```@docs
IntMult_Point_EllCurveZnZ(m::fmpz, P::Point_EllCurve{S}) where S <: Nemo.fmpz_mod
```

```@docs
rand_pair_EllCurve_Point(R::Oscar.MPolyRing_dec{S}, PP::Oscar.Geometry.ProjSpc{S}) where S <: Nemo.fmpz_mod
```

## Primality Proving

The Elliptic Curve Method (ECM) and the Elliptic Curve Primality Proving (ECPP) are methods involving elliptic curves. We introduce here these two functions and some related functions.  

### Methods

```@docs
ECM(n::fmpz; nbcurve::Int = 25000, multfact::fmpz = factorial(ZZ(10^4)))
```

```@docs
ECPP(n::fmpz)
```

```@docs
cornacchia_algorithm(d::fmpz, m::fmpz)
```

```@docs
Miller_Rabin_test(N::fmpz, k::Int64 = 20)
```

```@docs
Pollard_rho(N::fmpz, bound::Int = 50000)
```

```@docs
Pollard_p_1(N::fmpz, B::fmpz = ZZ(10)^5)
```
