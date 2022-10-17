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
DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{T}) where {T <: IntegerUnion}
ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}
```

### Addition, subtraction and scalar multiplication

Toric divisors can be added and subtracted via the usual `+` and `-`
operators. Moreover, multiplication by scalars from the left is supported
for scalars which are integers or of type `fmpz`.

### Special divisors

```@docs
trivial_divisor(v::AbstractNormalToricVariety)
anticanonical_divisor(v::AbstractNormalToricVariety)
canonical_divisor(v::AbstractNormalToricVariety)
```


## Properties of toric divisors

Equality of toric divisors can be tested via `==`.

To check if a toric divisor is trivial, one can invoke `is_trivial`.
This checks if all `coefficients` of the toric divisor in question
are zero. This must not be confused with a toric divisor being principal,
for which we support the following:
```@docs
is_principal(td::ToricDivisor)
```
Beyond this, we support the following properties of toric divisors:
```@docs
is_ample(td::ToricDivisor)
is_basepoint_free(td::ToricDivisor)
is_cartier(td::ToricDivisor)
is_effective(td::ToricDivisor)
is_integral(td::ToricDivisor)
is_nef(td::ToricDivisor)
is_prime(td::ToricDivisor)
is_q_cartier(td::ToricDivisor)
is_very_ample(td::ToricDivisor)
```


## Attributes

```@docs
coefficients(td::ToricDivisor)
polyhedron(td::ToricDivisor)
toric_variety(td::ToricDivisor)
```
