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
is_basepoint_free(td::ToricDivisor)
iscartier(td::ToricDivisor)
iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor)
isnef(td::ToricDivisor)
is_prime_divisor(td::ToricDivisor)
isprincipal(td::ToricDivisor)
is_q_cartier(td::ToricDivisor)
is_very_ample(td::ToricDivisor)
```

## Operations for toric divisors

```@docs
coefficients(td::ToricDivisor)
polyhedron(td::ToricDivisor)
toric_variety(td::ToricDivisor)
```
