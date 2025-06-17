```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Tropical semiring maps

## Introduction
In OSCAR, a `TropicalSemiringMap` is a map $\nu: K\to\mathbb{T}$ from a field $K$ to a tropical semiring $\mathbb{T}$ satisfying
1. finiteness: $\nu(a)=\pm\infty$ if and only if $a=0$,
2. multiplicativity: $\nu(a\cdot b)=\nu(a)+\nu(b)$,
3. superadditivity: $\nu(a\cdot b)\geq\min(\nu(a),\nu(b))$ (in the order defined in Section 2.7 of [Jos21](@cite)).

Most commonly, $\nu(a)=-\mathrm{val}(a)$ if $\mathbb{T}$ is the min-plus semiring, and $\nu(a)=+\mathrm{val}(a)$ if $\mathbb{T}$ is the max-plus semiring, for some valuation $\mathrm{val}:K^\ast\rightarrow\mathbb{R}$.  Essentially, $\nu$ captures a valuation on $K$ as well as a choice of min- or max-convention.  They are an optional input for most tropical functions over valued fields (the default being the trivial valuation and the min-convention).

## Constructor
Tropical semiring maps can be constructed as follows:
```@docs
tropical_semiring_map(K::Field, minOrMax::Union{typeof(min),typeof(max)}=min)
tropical_semiring_map(K::QQField, p::Union{RingElem,Integer,Rational}, minOrMax::Union{typeof(min),typeof(max)}=min)
tropical_semiring_map(Kt::Generic.RationalFunctionField, t::Generic.RationalFunctionFieldElem, minOrMax::Union{typeof(min),typeof(max)}=min)
```
