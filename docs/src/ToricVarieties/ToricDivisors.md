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

### Special constructors

Addition of toric divisors `td1` and `td2` (on the same toric variety) and
scalar multiplication with `c` (it can be either valued in `Int64` or `fmpz`)
is supported via `c * td1 + td2`. One can subtract them via `td1 - td2`.


### Equality

Equality of two toric divisors `td1` and `td2` (on the same toric variety)
is achieved by checking if their coefficients are identical.
This is implemented via `td1 == td2`.


## Properties of toric divisors

To check if a toric divisor `td` is trivial, one can invoke `istrivial(td)`.
Internally, this executes the following method:
```@docs
isprincipal(td::ToricDivisor)
```
Beyond this, we support the following properties of toric divisors:
```@docs
isample(td::ToricDivisor)
is_basepoint_free(td::ToricDivisor)
iscartier(td::ToricDivisor)
iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor)
isnef(td::ToricDivisor)
isprime(td::ToricDivisor)
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
