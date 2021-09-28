```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["tg_divisors.md"]
```

# Toric Divisors


## Introduction

Toric divisors are those divisors that are invariant under the torus action.
They are formal sums of the codimension one orbits, and these in turn
correspond to the rays of the underlying fan.


## Construction

```@docs
ToricDivisor(coeffs::AbstractVector, X::NormalToricVarietyType)
```

## Auxiliary functions
```@docs
isample(td::ToricDivisor)
isbasepoint_free(td::ToricDivisor)
iscartier(td::ToricDivisor)
toric_iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor)
isnef(td::ToricDivisor)
isprincipal(td::ToricDivisor)
isq_cartier(td::ToricDivisor)
issemiample(td::ToricDivisor)
isvery_ample(td::ToricDivisor)
polyhedron_of_divisor(td::ToricDivisor)
```
