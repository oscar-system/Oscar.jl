```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Abstract Variety Maps

## Constructors

```@docs
map(X::AbstractVariety, Y::AbstractVariety, fˣ::Vector, fₓ = nothing; inclusion::Bool = false, symbol::String = "x")
```

```@docs
identity_map(X::AbstractVariety)
```

## Underlying Data of an Abstract Variety Map

An abstract variety map is made up from (a selection of) the data discussed here:

```@docs
domain(f:: AbstractVarietyMap)
```

```@docs
codomain(f:: AbstractVarietyMap)
```

```@docs
dim(f::AbstractVarietyMap)
```

```@docs
pullback(f::AbstractVarietyMap, y::MPolyDecRingElem)
```

```@docs
pushforward(f::AbstractVarietyMap, x::MPolyDecRingElem)
```

```@docs
tangent_bundle(f::AbstractVarietyMap)
```

## Further Data Associated to an Abstract Variety Map

```@docs
cotangent_bundle(f::AbstractVarietyMap)
```

```@docs
todd_class(f::AbstractVarietyMap)
```

## Operations on Abstract Variety Maps

```@docs
compose(f::AbstractVarietyMap, g::AbstractVarietyMap)
```

