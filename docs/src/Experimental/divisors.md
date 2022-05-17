```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["divisors.md"]
```

# Divisors

In order to consider divisors on curves, we restrict our attention to smooth and irreducible curves.

Let ``C`` be an affine or projective plane curve defined by an irreducible equation ``F``. Then any polynomial function ``G`` which is not divisible by ``F`` will vanish on ``C`` only at finitely many points. A way to encode these points together with their intersection multiplicities is to consider a divisor. A divisor on a curve is a formal finite sum of points of the curve with integer coefficients. A natural operation of addition can be defined on the set of divisors of a curve, which turns it into an Abelian group.


## Constructors

Divisors on curves are here introduced as a dictionary associating a point on the curve to its multiplicity.

```@docs
AffineCurveDivisor
```

```@docs
ProjCurveDivisor
```

To define the divisor ``0`` of the group of divisors, one uses the following method:

```@docs
curve_zero_divisor(C::ProjPlaneCurve{S}) where S <: FieldElem
curve_zero_divisor(C::AffinePlaneCurve{S}) where S <: FieldElem
```

## Methods

The following functions on divisors of curves are implemented.

```@docs
curve(D::Oscar.PlaneCurveModule.CurveDivisor)
degree(D::Oscar.PlaneCurveModule.CurveDivisor)
is_effective(D::Oscar.PlaneCurveModule.CurveDivisor)
is_linearly_equivalent(D::ProjCurveDivisor, E::ProjCurveDivisor)
is_principal(D::ProjCurveDivisor{S}) where S <: FieldElem
principal_divisor(D::ProjCurveDivisor{S}) where S <: FieldElem
global_sections(D::ProjCurveDivisor)
dimension_global_sections(D::ProjCurveDivisor)
```


In addition, the multiplicity of a polynomial or a fraction at a given point can be computed:

```@docs
multiplicity(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}, P::Point{S}) where {S <: FieldElem, T <: MPolyElem{S}}
multiplicity(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}, P::Point{S}) where S <: FieldElem
```


The divisor of a polynomial or a fraction along a plane curve can be computed, but will give only the points which belong to the base field.


```@docs
divisor(C::AffinePlaneCurve{S}, F::Oscar.MPolyElem{S}) where S <: FieldElem
divisor(C::AffinePlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T}) where {S <: FieldElem, T <: MPolyElem{S}}
divisor([PP::Oscar.Geometry.ProjSpc{S}], C::ProjPlaneCurve{S}, F::Oscar.MPolyElem_dec{S}) where S <: FieldElem
divisor(PP::Oscar.Geometry.ProjSpc{S}, C::ProjPlaneCurve{S}, phi::AbstractAlgebra.Generic.Frac{T})  where {S <: FieldElem, T <: Oscar.MPolyElem_dec{S}}
```
