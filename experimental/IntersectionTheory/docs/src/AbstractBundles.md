```@meta
CurrentModule = Oscar
```

# Abstract Bundles

## Constructors

```@docs
abstract_bundle(X::AbstractVariety, ch::Union{MPolyDecRingElem, MPolyQuoRingElem})
```

## Underlying Data of an Abstract Bundle

```@docs
parent(F::AbstractBundle)
```

```@docs
rank(F::AbstractBundle)
```

```@docs
chern_character(F::AbstractBundle)
```

```@docs
total_chern_class(F::AbstractBundle)
```

## Further Data Associated to an Abstract Bundle

```@docs
chern_class(F::AbstractBundle, k::Int)
```

```@docs
top_chern_class(F::AbstractBundle)
```

```@docs
segre_class(F::AbstractBundle)
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


## Operations on Abstract Bundles

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
