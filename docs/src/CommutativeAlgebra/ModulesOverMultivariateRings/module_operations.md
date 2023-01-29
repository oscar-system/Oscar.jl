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
kernel(a::ModuleFPHom)
```

### Image

```@docs
image(a::ModuleFPHom)
```

### Cokernel

```@docs
cokernel(a::ModuleFPHom)
```

## Direct Sums and Products

```@docs
direct_sum(M::ModuleFP{T}...; task::Symbol = :sum) where T
```

```@docs
direct_product(M::ModuleFP{T}...; task::Symbol = :prod) where T
```
