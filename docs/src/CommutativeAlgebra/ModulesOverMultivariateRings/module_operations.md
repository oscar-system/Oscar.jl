```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
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
direct_sum(M::ModuleFP{T}, Ms::ModuleFP{T}...; task::Symbol = :sum) where T
```

```@docs
direct_product(M::ModuleFP{T}, Ms::ModuleFP{T}...; task::Symbol = :prod) where T
```

## Truncation

```@docs
truncate(M::ModuleFP, g::FinGenAbGroupElem, task::Symbol=:with_morphism)
```

## Twists

In the graded case, we have:

```@docs
twist(M::ModuleFP{T}, g::FinGenAbGroupElem) where {T<:MPolyDecRingElem}
```
