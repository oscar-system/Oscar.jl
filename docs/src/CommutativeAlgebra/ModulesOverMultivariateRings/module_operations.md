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

## Presentations

```@docs
presentation(M::ModuleFP)
```

## Syzygies and Free Resolutions

```@docs
free_resolution(M::SubQuo; ordering::ModuleOrdering = default_ordering(M),
    length::Int=0, algorithm::Symbol=:fres)
```

## Homology

```@docs
homology(C::ChainComplex{<:ModuleFP})
```

```@docs
homology(C::ChainComplex{<:ModuleFP}, i::Int)
```

## Hom and Ext

```@docs
hom(M::ModuleFP, N::ModuleFP, alg::Symbol=:maps)
```

```@docs
element_to_homomorphism(f::ModuleFPElem)
```

```@docs
homomorphism_to_element(H::ModuleFP, phi::ModuleFPHom)
```

```@docs
ext(M::ModuleFP, N::ModuleFP, i::Int)
```

## Tensorproduct and Tor

```@docs
tensor_product(G::ModuleFP...; task::Symbol = :none)
```

```@docs
tor(M::ModuleFP, N::ModuleFP, i::Int)
```

