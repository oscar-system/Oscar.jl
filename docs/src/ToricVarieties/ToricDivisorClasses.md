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

```@docs
Base.:+(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
Base.:-(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
Base.:*(c::fmpz, td::ToricDivisorClass)
```

### Equality

```@docs
Base.:(==)(tdc1::ToricDivisorClass, tdc2::ToricDivisorClass)
```


## Properties of toric divisor classes

```@docs
istrivial(tdc::ToricDivisorClass)
```


## Operations for toric divisor classes

```@docs
divisor_class(tdc::ToricDivisorClass)
toric_variety(tdc::ToricDivisorClass)
toric_divisor(tdc::ToricDivisorClass)
```
