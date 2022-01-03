```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["ToricLineBundles.md"]
```


# Toric Line Bundles


## Constructors

```@docs
ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{fmpz})
ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{Int})
```


## Properties

```@docs
istrivial(l::ToricLineBundle)
is_basepoint_free(l::ToricLineBundle)
isample(l::ToricLineBundle)
is_very_ample(l::ToricLineBundle)
```


## Attributes

```@docs
degree(l::ToricLineBundle)
divisor_class(l::ToricLineBundle)
toric_divisor(l::ToricLineBundle)
variety(l::ToricLineBundle)
```
