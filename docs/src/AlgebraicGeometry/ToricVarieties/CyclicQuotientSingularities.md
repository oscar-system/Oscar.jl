```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# Cyclic Quotient Singularities

## Introduction

Cyclic quotient singularities are quotients of $\mathbb{C}^2$ by the action of
$\mathbb{Z}/n\mathbb{Z}$ acting via 
$$\left(\begin{array}{cc}\xi & 0\\0 & \xi^q\end{array}\right)$$,
where $\xi$ is a $n$-th root of unity, and $q$ and $n$ are integers, such that $q$ is coprime with $n$, and $0<q<n$.

For the notation we rely on [Chr91](@cite) and [Ste91](@cite).

!!! warning
    Note that [Chr91](@cite) and [Ste91](@cite) use Hirzebruch-Jung continued fraction, which differ from the
    commonly known continued fraction from literature and used in the rest of OSCAR.


## Constructors

```@docs
cyclic_quotient_singularity(n::T, q::T) where {T <: IntegerUnion}
```


## Attributes

```@docs
continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
```


## Auxiliary Methods

```@docs
continued_fraction_hirzebruch_jung_to_rational(v::Vector{ZZRingElem})
rational_to_continued_fraction_hirzebruch_jung(r::QQFieldElem)
```
