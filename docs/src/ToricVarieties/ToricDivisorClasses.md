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
ToricDivisorClass(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}
```

### Addition, subtraction and scalar multiplication

Toric divisor classes can be added and subtracted via the usual `+` and `-`
operators. Moreover, multiplication by scalars from the left is supported
for scalars which are integers or of type `fmpz`.

### Special divisor classes

```@docs
trivial_divisor_class(v::AbstractNormalToricVariety)
anticanonical_divisor_class(v::AbstractNormalToricVariety)
canonical_divisor_class(v::AbstractNormalToricVariety)
```


## Properties

Equality of toric divisor classes can be tested via `==`.

To check if a toric divisor class is trivial, one can invoke `is_trivial`.


## Attributes

```@docs
divisor_class(tdc::ToricDivisorClass)
toric_variety(tdc::ToricDivisorClass)
toric_divisor(tdc::ToricDivisorClass)
```
