```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Rational Points on Projective Schemes

```@docs
AbsProjectiveRationalPoint
ProjectiveRationalPoint
rational_point(X::AbsProjectiveScheme, coordinates; check::Bool=true)
coordinates(P::ProjectiveRationalPoint)
ideal(P::AbsProjectiveRationalPoint)
scheme(P::ProjectiveRationalPoint)
normalize!(a::AbsProjectiveRationalPoint{<:FieldElem})
normalize!(a::AbsProjectiveRationalPoint{ZZRingElem})
```
