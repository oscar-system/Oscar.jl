```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Projective

```@docs
AbsProjectiveRationalPoint
ProjectiveRationalPoint
coordinates(P::ProjectiveRationalPoint)
ideal(P::ProjectiveRationalPoint)
scheme(P::ProjectiveRationalPoint)
normalize!(a::AbsProjectiveRationalPoint{<:FieldElem})
normalize!(a::AbsProjectiveRationalPoint{ZZRingElem})
```
