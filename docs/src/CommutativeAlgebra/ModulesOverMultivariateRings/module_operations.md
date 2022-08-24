```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["module_operations.md"]
```

# Operations on Modules

## Subquotients Related to Homomorphisms

### Kernel

```@docs
kernel(a::FreeModuleHom)
```

```@docs
kernel(a::SubQuoHom)
```


### Image

```@docs
image(a::FreeModuleHom)
```

```@docs
image(a::SubQuoHom)
```

### Cokernel

```@docs
cokernel(a::FreeModuleHom) 
```

## Chain Complexes

### Constructors

### Data Associated to Chain Complexes

### Operations on Chain Complexes

### Tests on Chain Complexes

### Homology


## Direct Sums and Products

```@docs
direct_sum(M::ModuleFP{T}...; task::Symbol = :none) where T
```

```@docs
direct_product(M::ModuleFP{T}...; task::Symbol = :none) where T
```

```@docs
hom_product(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleMap})
```

## Presentations


## Syzygies and Free Resolutions

```@docs
free_resolution(M::SubQuo; ordering::ModuleOrdering = default_ordering(M),
    length::Int=0, algorithm::Symbol=:fres)
```


## Hom and Ext



## Tensorproduct and Tor

```@docs
hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{ <: ModuleMap})
```
