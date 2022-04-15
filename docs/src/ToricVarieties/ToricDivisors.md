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

### General constructors

```@docs
DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{Int})
ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{Int})
```

### Special constructors

```@docs
Base.:+(td1::ToricDivisor, td2::ToricDivisor)
Base.:-(td1::ToricDivisor, td2::ToricDivisor)
Base.:*(c::fmpz, td::ToricDivisor)
```

### Equality

```@docs
Base.:(==)(td1::ToricDivisor, td2::ToricDivisor)
```


## Properties of toric divisors

```@docs
isample(td::ToricDivisor)
is_basepoint_free(td::ToricDivisor)
iscartier(td::ToricDivisor)
iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor)
isnef(td::ToricDivisor)
isprime(td::ToricDivisor)
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

## Special divisors

```@docs
trivial_divisor(v::AbstractNormalToricVariety)
anticanonical_divisor(v::AbstractNormalToricVariety)
canonical_divisor(v::AbstractNormalToricVariety)
```
