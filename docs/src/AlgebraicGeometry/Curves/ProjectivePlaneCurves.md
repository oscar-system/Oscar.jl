```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```
# Projective Plane Curves
```@docs
ProjectivePlaneCurve
```

Projective plane curves are modeled in Oscar as projective
algebraic sets. See `AbsProjectiveAlgebraicSet`(@ref).
In addition to the methods for algebraic sets and curves
the following methods special to plane curves are available.

```@docs
defining_equation(C::ProjectivePlaneCurve{S,MPolyQuoRing{T}}) where {S,T}
degree(C::ProjectivePlaneCurve)
common_components(C::S, D::S) where {S<:ProjectivePlaneCurve}
multiplicity(C::ProjectivePlaneCurve, P::AbsProjectiveRationalPoint)
tangent_lines(C::ProjectivePlaneCurve, P::AbsProjectiveRationalPoint)
intersection_multiplicity(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve
is_transverse_intersection(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve
```

