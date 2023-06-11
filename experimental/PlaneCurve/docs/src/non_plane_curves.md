```@meta
CurrentModule = Oscar
```


# Projective Curves

We consider projective curves in projective spaces of arbitrary dimension.

## Constructors

We define a projective curve by an ideal of homogeneous polynomials.

```@docs
ProjCurve
```

## General functions for curves

```@docs
defining_ideal(C::ProjCurve)
in(P::Oscar.Geometry.ProjSpcElem, C::ProjCurve)
curve_components(C::ProjCurve)
is_irreducible(C::ProjCurve)
reduction(C::ProjCurve)
jacobi_ideal(C::ProjCurve)
invert_birational_map(phi::Vector{T}, C::ProjCurve) where {T <: MPolyRingElem}
```
