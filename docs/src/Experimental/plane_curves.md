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
Pages = ["plane_curves.md"]
```


# Affine and Projective Plane Curves

We consider two kinds of plane curves: affine and projective. An affine plane
curve is defined by a polynomial in two variables, whereas a projective plane
curve is defined by a homogeneous polynomial belonging to a graded polynomial
ring in three variables.

## Affine Plane Curves

An affine plane curve is defined as the class of a two-variables polynomial
$F$ over a field $K$, modulo the equivalence relation $F \sim G \iff
\exists \lambda \in K\backslash \{0\}, F = \lambda \cdot G$.

### Example

```@docs
AffinePlaneCurve
```

## Projective Plane Curves

Similarly, a projective plane curve is defined as the class of a
three-variables homogeneous polynomial $F$ over a field $K$, modulo the
equivalence relation $F\sim G \iff \exists \lambda \in K\backslash \{0\}, F =
\lambda \cdot G$. The defining equation is supposed to belong to a graded
ring.

```@docs
ProjPlaneCurve
```

A particular kind of projective curves is the case of elliptic curves, see the
corresponding section for more information. The types `ProjPlaneCurve` and
`ProjEllipticCurve` are subtypes of the abstract type `ProjectivePlaneCurve`.
In addition, the types `AffinePlaneCurve` and `ProjectivePlaneCurve` are
subtypes of the abstract type `PlaneCurve`.

## Points

When considering curves, it is natural to have a look at points on the curve.
We describe in this section how to deal with points, both in the affine and
projective settings.

### Point in the affine space

A point in the affine space can be defined as follows:
```@docs
Point
```
We consider also the following function for points.

```@docs
ideal_point(R::MPolyRing{S}, P::Point{S}) where S <: FieldElem
```

The following function checks if a given point is on a curve:

```@docs
in(P::Point{S}, C::AffinePlaneCurve{S}) where S <: FieldElem
```

### Point in the projective space

In order to define a point in the projective plane, one needs first to define
the projective plane as follows, where `K` is the base ring:

#### Example
```jldoctest
julia> K = QQ
Rational Field

julia> PP = proj_space(K, 2)
(Projective space of dim 2 over Rational Field
, MPolyElem_dec{fmpq, fmpq_mpoly}[x[0], x[1], x[2]])

```

Then, one can define a projective point as follows:

#### Example
```jldoctest
julia> P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(2), QQ(-5)])
ERROR: UndefVarError: PP not defined
Stacktrace:
 [1] top-level scope
   @ none:1

```

The following function checks if a given point is on a curve:

```@docs
in(P::Oscar.Geometry.ProjSpcElem{S}, C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
```

## General functions for curves

```@docs
degree(C::Oscar.PlaneCurveModule.PlaneCurve)
ring(C::Oscar.PlaneCurveModule.PlaneCurve)
curve_components(C::Oscar.PlaneCurveModule.PlaneCurve{S}) where S <: FieldElem
reduction(C::AffinePlaneCurve{S}) where S <: FieldElem
is_irreducible(C::Oscar.PlaneCurveModule.PlaneCurve{S}) where S <: FieldElem
is_reduced(C::Oscar.PlaneCurveModule.PlaneCurve{S}) where S <: FieldElem
union(C::T, D::T) where T <: Oscar.PlaneCurveModule.PlaneCurve
arithmetic_genus(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve)
arithmetic_genus(C::AffinePlaneCurve)
geometric_genus(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
geometric_genus(C::AffinePlaneCurve)
```

## Smoothness, tangents and singularity related functions

```@docs
jacobi_ideal(C::Oscar.PlaneCurveModule.PlaneCurve)
is_smooth(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
is_smooth(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
tangent(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
tangent(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
curve_singular_locus(C::AffinePlaneCurve)
curve_singular_locus([PP::Oscar.Geometry.ProjSpc{S}], C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
multiplicity(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
multiplicity(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
tangent_lines(C::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
tangent_lines(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
is_smooth_curve(C::AffinePlaneCurve)
is_smooth_curve(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve)
```

## Intersection of curves

```@docs
common_components(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
common_components(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
curve_intersect(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}) where S <: FieldElem
curve_intersect([PP::Oscar.Geometry.ProjSpc{S}], C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}) where S <: FieldElem
intersection_multiplicity(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S <: FieldElem
intersection_multiplicity(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S <: FieldElem
aretransverse(C::AffinePlaneCurve{S}, D::AffinePlaneCurve{S}, P::Point{S}) where S<:FieldElem
aretransverse(C::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, D::Oscar.PlaneCurveModule.ProjectivePlaneCurve{S}, P::Oscar.Geometry.ProjSpcElem{S}) where S<:FieldElem
```

