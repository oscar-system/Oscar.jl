```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract variety maps

## Types

The OSCAR type for abstract variety maps is `AbstractVarietyMap`.

## Constructors

```@docs
map(X::AbstractVariety, Y::AbstractVariety, fˣ::Vector, fₓ = nothing; inclusion::Bool = false, symbol::String = "x")
```

```@docs
identity_map(X::AbstractVariety)
```

## Underlying data of an abstract variety map

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
pullback(f::AbstractVarietyMap, y::MPolyDecRingOrQuoElem)
```

```@docs
pushforward(f::AbstractVarietyMap, x::MPolyDecRingOrQuoElem)
```

```@docs
tangent_bundle(f::AbstractVarietyMap)
```

## Further data associated to an abstract variety map

```@docs
cotangent_bundle(f::AbstractVarietyMap)
```

```@docs
todd_class(f::AbstractVarietyMap)
```

## Operations on abstract variety maps

```@docs
compose(f::AbstractVarietyMap, g::AbstractVarietyMap)
```

