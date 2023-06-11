```@meta
CurrentModule = Oscar
```


# Introduction

This page lists the OSCAR features for tropical geometry, which are still at the very beginning of their development. For notation we refer to [MS15](@cite) and [Jos21](@cite).


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Lars Kastner](https://lkastner.github.io/),
* [Marta Panizzut](https://martapanizzut.github.io/),
* [Yue Ren](https://www.yueren.de/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).


## Tropical semirings

```@docs
TropicalSemiring
TropicalSemiringMap
tropical_polynomial
```


## Tropical curves

(no functions yet, under active development)

```@docs
TropicalCurve
```

in work: chip firing, jacobians.


## Tropical hypersurfaces

```@docs
TropicalHypersurface
dual_subdivision(TH::TropicalHypersurface{M, EMB}) where {M, EMB}
```

in work: minimal tropical polynomials of a hypersurface.


## Tropical linear spaces

```@docs
TropicalLinearSpace
```


## Tropical varieties

(no functions yet, under active development)

in work: `tropical_variety(I::MPolyIdeal,val::TropicalSemiringMap)`


## Groebner bases and Groebner polyhedra

```@docs
groebner_basis(I::MPolyIdeal,val::TropicalSemiringMap,w::Vector{<: Union{Int,Rational{Int}} }; complete_reduction::Bool=false, return_lead::Bool=false)
```

in work: `groebner_polyhedron(I::MPolyIdeal,val::TropicalSemiringMap,w::Vector{<: Union{Int,Rational{Int}})`


## Intersections and stable intersections


```@docs
intersect(T1::TropicalVarietySupertype{M, EMB}, T2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
stably_intersect(T1::TropicalVarietySupertype{M, EMB}, T2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
```
