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



## Construction

```@docs
ToricDivisor(coeffs::AbstractVector, X::NormalToricVarietyType)
```

## Auxiliary functions
```@docs
isample(td::ToricDivisor)
isbasepoint_free(td::ToricDivisor)
iscartier(td::ToricDivisor)
iseffective(td::ToricDivisor)
isintegral(td::ToricDivisor)
isnef(td::ToricDivisor)
isprincipal(td::ToricDivisor)
isq_cartier(td::ToricDivisor)
issemiample(td::ToricDivisor)
isvery_ample(td::ToricDivisor)
```
