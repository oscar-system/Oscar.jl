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

## Attributes

```@docs
divisor_class(l::ToricLineBundle)
variety(l::ToricLineBundle)
```
