```@meta
CurrentModule = Oscar
```

# Abstract Variety Morphisms

## Constructors

```@docs
hom(X::AbstractVariety, Y::AbstractVariety, fˣ::Vector, fₓ = nothing; inclusion::Bool = false, symbol::String = "x")
```

## Underlying Data of an Abstract Variety

domain::V1
  codomain::V2
  dim::Int
  pullback::AffAlgHom
  pushforward::FunctionalMap
  O1::MPolyDecRingOrQuoElem
  T::AbstractBundle{V1}

## Further Data Associated to an Abstract Variety


## Operations on Abstract Varietiy Maps
