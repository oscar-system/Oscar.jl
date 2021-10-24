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
ToricDivisor( v::AbstractNormalToricVariety, coeffs::Vector{Int} )
DivisorOfCharacter( v::AbstractNormalToricVariety, character::Vector{Int} )
```

## Properties of toric divisors

```@docs
isample
isbasepoint_free
iscartier
iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor) 
isnef
isprime_divisor
isprincipal
isq_cartier
isvery_ample
```

## Operations for toric divisors

```@docs
polyhedron(td::ToricDivisor)
coefficients(td::ToricDivisor)
```
