```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Tropical semirings, matrices, and polynomials

## Introduction
In OSCAR, the tropical semiring is either
- the min-plus semiring $(\mathbb{Q}\cup\{+\infty\},\oplus,\odot)$ with $a\oplus b=\min(a,b)$ and $a\odot b=a+b$,
- the max-plus semiring $(\mathbb{Q}\cup\{-\infty\},\oplus,\odot)$ with $a\oplus b=\max(a,b)$ and $a\odot b=a+b$.

Whereas tropical semirings in [MS15](@cite) and [Jos21](@cite) are extensions of the real numbers, tropical semirings in OSCAR are an extension of the rational numbers to avoid precision issues.

## Constructor
Objects of type `TropicalSemiring`, as well as matrices and polynomials thereover, can be constructed as follows:
```@docs
tropical_semiring()
```

## Properties
Objects of type `TropicalSemiring` have the following properties:
```@docs
convention(T::TropicalSemiring{typeof(min)})
```

## Related functions
Other functions related to `TropicalSemiring`, matrices, and polynomials thereover include:
```@docs
det(A::AbstractAlgebra.Generic.MatSpaceElem{<: Oscar.TropicalSemiringElem})
roots(f::PolyRingElem{TropicalSemiringElem{minOrMax}}) where {minOrMax <: Union{typeof(min),typeof(max)}}
tropical_polynomial(f::MPolyRingElem, nu::TropicalSemiringMap)
```
