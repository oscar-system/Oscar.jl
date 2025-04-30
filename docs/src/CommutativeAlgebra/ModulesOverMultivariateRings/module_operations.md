```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# Operations on Modules

## Subquotients Related to Homomorphisms

### Kernel

```@docs
kernel(a::SparseFPModuleHom)
```

### Image

```@docs
image(a::SparseFPModuleHom)
```

### Cokernel

```@docs
cokernel(a::SparseFPModuleHom)
```

## Direct Sums and Products

```@docs
direct_sum(M::SparseFPModule{T}, Ms::SparseFPModule{T}...; task::Symbol = :sum) where T
```

```@docs
direct_product(M::SparseFPModule{T}, Ms::SparseFPModule{T}...; task::Symbol = :prod) where T
```

## Truncation

```@docs
truncate(M::SparseFPModule, g::FinGenAbGroupElem, task::Symbol=:with_morphism)
```

## Twists

In the graded case, we have:

```@docs
twist(M::SparseFPModule{T}, g::FinGenAbGroupElem) where {T<:MPolyDecRingElem}
```
