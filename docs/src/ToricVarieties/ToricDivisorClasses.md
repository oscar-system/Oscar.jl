```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["ToricDivisorClasses.md"]
```


# Toric Divisor Classes

## Introduction

Toric divisor classes are equivalence classes of Weil divisors modulo linear equivalence.


## Constructors

### General constructors

```@docs
ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{fmpz})
ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{Int})
```

### Special constructors

Addition of toric divisor classes `tdc1` and `tdc2` (on the same toric variety) and
scalar multiplication with `c` (it can be either valued in `Int64` or `fmpz`)
is supported via `c * tdc1 + tdc2`. One can subtract them via `tdc1 - tdc2`.


### Equality

Equality of two toric divisor classes `tdc1` and `tdc2` (on the same toric variety)
is achieved by checking if their difference is a trivial class, i.e. the divisor class of
a principal toric divisor. This is implemented via `tdc1 == tdc2`.


## Properties of toric divisor classes

To check if a toric divisor class `tdc` is trivial, one can invoke `is_trivial(tdc)`.


## Operations for toric divisor classes

```@docs
divisor_class(tdc::ToricDivisorClass)
toric_variety(tdc::ToricDivisorClass)
toric_divisor(tdc::ToricDivisorClass)
```

## Special divisor classes

```@docs
trivial_divisor_class(v::AbstractNormalToricVariety)
anticanonical_divisor_class(v::AbstractNormalToricVariety)
canonical_divisor_class(v::AbstractNormalToricVariety)
```
