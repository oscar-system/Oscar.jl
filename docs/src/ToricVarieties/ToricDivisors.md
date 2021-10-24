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
ToricDivisor
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
isprincipal
isq_cartier
isvery_ample
polyhedron(td::ToricDivisor)
```

