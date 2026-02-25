```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract bundles

## Types

The OSCAR type for abstract vector bundles is `AbstractBundle`.

## Constructors

```@docs
abstract_bundle(X::AbstractVariety, ch::Union{MPolyDecRingElem, MPolyQuoRingElem})
```

## Underlying data of an abstract bundle

An abstract bundle is made up from (a selection of) the data discussed here:

```@docs
parent(F::AbstractBundle)
```

```@docs
chern_character(F::AbstractBundle)
```

```@docs
rank(F::AbstractBundle)
```

```@docs
total_chern_class(F::AbstractBundle)
```

## Further data associated to an abstract bundle

```@docs
chern_class(F::AbstractBundle, k::Int)
```

```@docs
top_chern_class(F::AbstractBundle)
```

```@docs
total_segre_class(F::AbstractBundle)
```

```@docs
segre_class(F::AbstractBundle, k::Int)
```

```@docs
 todd_class(F::AbstractBundle)
```

```@docs
total_pontryagin_class(F::AbstractBundle)
```

```@docs
pontryagin_class(F::AbstractBundle, k::Int)
```

```@docs
euler_characteristic(F::AbstractBundle)
```

```@docs
hilbert_polynomial(F::AbstractBundle)
```

## Operations on abstract bundles

```@docs
-(F::AbstractBundle)
```

```@docs
dual(F::AbstractBundle)
```

```@docs
det(F::AbstractBundle)
```

```@docs
exterior_power(F::AbstractBundle, k::Int)
```

```@docs
symmetric_power(F::AbstractBundle, k::Int)
```

```@docs
pullback(f::AbstractVarietyMap, F::AbstractBundle)
```

```@docs
pushforward(f::AbstractVarietyMap, F::AbstractBundle)
```

## Tests on abstract bundles

```@docs
==(F::AbstractBundle, G::AbstractBundle)
```
