```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["GeneralSchemes.md"]
```

# General schemes

Arbitrary schemes over a commutative base ring ``\mathbb k`` with unit 
are instances of the abstract type
```@docs
Scheme{BaseRingType<:Ring}
```
Morphisms of schemes shall be derived from the abstract type
```@docs
    SchemeMor{DomainType, CodomainType, MorphismType, BaseMorType}
```
