```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Rational Points on Projective Schemes

```@docs
AbsProjectiveRationalPoint
ProjectiveRationalPoint
coordinates(P::ProjectiveRationalPoint)
ideal(P::AbsProjectiveRationalPoint)
scheme(P::ProjectiveRationalPoint)
normalize!(a::AbsProjectiveRationalPoint{<:FieldElem})
normalize!(a::AbsProjectiveRationalPoint{ZZRingElem})
```
