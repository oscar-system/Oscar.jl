```@meta
CurrentModule = Oscar
```

# Affine plane curves

```@docs
AffinePlaneCurve
defining_equation(C::AffinePlaneCurve{S, MPolyQuoRing{E}}) where {S, E}
common_components(C::S, D::S) where {S<:AffinePlaneCurve}
multiplicity(C::AffinePlaneCurve, P::AbsAffineRationalPoint)
tangent_lines(C::AffinePlaneCurve, P::AbsAffineRationalPoint)
intersection_multiplicity(C::AffinePlaneCurve, D::AffinePlaneCurve, P::AbsAffineRationalPoint)
is_transverse_intersection(C::AffinePlaneCurve, D::AffinePlaneCurve, P::AbsAffineRationalPoint)
projective_closure(C::AffinePlaneCurve)
```
