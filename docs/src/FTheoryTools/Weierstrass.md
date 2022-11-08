```@meta
CurrentModule = FTheoryTools
```

```@contents
Pages = ["Weierstrass.md"]
```

# Global Weierstrass models

## Constructors

We support the following constructors:
```@docs
GlobalWeierstrassModel(polys::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}})
GenericGlobalWeierstrassModelOverToricSpace(v::Oscar.AbstractNormalToricVariety)
GenericGlobalWeierstrassModelOverProjectiveSpace(n::Int)
```

## Attributes

```@docs
poly_f(w::GlobalWeierstrassModel)
poly_g(w::GlobalWeierstrassModel)
toric_base_space(w::GlobalWeierstrassModel)
```
