```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["ToricDivisors.md"]
```


# Toric Divisors

## Introduction

Toric divisors are those divisors that are invariant under the torus action.
They are formal sums of the codimension one orbits, and these in turn
correspond to the rays of the underlying fan.


## Constructors

```@docs
DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{Int})
ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{Int})
```

## Properties of toric divisors

```@docs
isample(td::ToricDivisor)
isbasepoint_free(td::ToricDivisor)
iscartier(td::ToricDivisor)
iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor)
isnef(td::ToricDivisor)
isprime_divisor(td::ToricDivisor)
isprincipal(td::ToricDivisor)
isq_cartier(td::ToricDivisor)
isvery_ample(td::ToricDivisor)
```

## Operations for toric divisors

```@docs
coefficients(td::ToricDivisor)
polyhedron(td::ToricDivisor)
toricvariety(td::ToricDivisor)
```
