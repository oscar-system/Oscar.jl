```@meta
CurrentModule = Oscar
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

## Change of base
```@docs
base_change(phi::Any, X::Scheme)
base_change(phi::Any, f::SchemeMor;
    domain_map::AbsSchemeMor, codomain_map::AbsSchemeMor
    )
```
   
