```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["standard_constructions.md"]
```

# Standard Constructions in Algebraic Geometry

This page documents a collection of standard constructions in algebraic geometry
available in Oscar.

## Grassman Plücker Ideal
```@docs
grassman_pluecker_ideal([ring::MPolyRing,] subspace_dimension::Int, ambient_dimension::Int)
```

