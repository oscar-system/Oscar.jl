```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Rational Points on Affine Schemes

```@docs
AbsAffineRationalPoint
AffineRationalPoint
rational_point(X::AbsAffineScheme, coordinates; check::Bool=true)
coordinates(p::AffineRationalPoint)
ideal(P::AbsAffineRationalPoint)
scheme(P::AbsAffineRationalPoint)
closed_embedding(P::AbsAffineRationalPoint)
is_smooth(P::AbsAffineRationalPoint)
tangent_space(P::AbsAffineRationalPoint{<:FieldElem})
```

Some experimental methods are available too.
Note that their interface is likely to change in the future.

```@docs
is_du_val_singularity(P::AbsAffineRationalPoint{<:FieldElem})
decide_du_val_singularity(P::AbsAffineRationalPoint{<:FieldElem,<:Any})
```
