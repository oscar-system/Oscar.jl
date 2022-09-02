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
kernel(a::ModuleMap)
```

### Image

```@docs
image(a::ModuleMap)
```

### Cokernel

```@docs
cokernel(a::ModuleMap)
```

## Chain Complexes

### Constructors

### Data Associated to Chain Complexes

### Operations on Chain Complexes

### Tests on Chain Complexes

### Maps of Chain Complexes

#### Types

#### Constructors

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
homomorphism(f::ModuleFPElem)
```

```@docs
module_elem(H::ModuleFP, phi::ModuleMap)
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

