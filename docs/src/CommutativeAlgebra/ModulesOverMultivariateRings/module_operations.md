```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Operations on Modules

## Subquotients Related to Homomorphisms

### Kernel

```@docs
kernel(a::OFPModuleHom)
```

### Image

```@docs
image(a::OFPModuleHom)
```

### Cokernel

```@docs
cokernel(a::OFPModuleHom)
```

## Direct Sums and Products

```@docs
direct_sum(M::OFPModule{T}, Ms::OFPModule{T}...; task::Symbol = :sum) where T
```

```@docs
direct_product(M::OFPModule{T}, Ms::OFPModule{T}...; task::Symbol = :prod) where T
```

## Annihilator

```@docs
annihilator(N::OFPModule{T}) where T
```

## Truncation

```@docs
truncate(M::OFPModule, g::FinGenAbGroupElem, task::Symbol=:with_morphism)
```

## Twists

In the graded case, we have:

```@docs
twist(M::OFPModule{T}, g::FinGenAbGroupElem) where {T<:MPolyDecRingElem}
```
