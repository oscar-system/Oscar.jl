```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["homological_algebra.md"]
```

# Homological Algebra

Some OSCAR functions which are fundamental to homological algebra such as the `kernel` function
and basic functions for handling chain and cochain complexes have already been discussed
in previous sections. Building on these functions, we now introduce further OSCAR functionality
supporting computations in homological algebra.

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

## Flatness

```@julia
fitting_ideal(M::ModuleFP, i::Int)
```

```@julia
is_flat(M::ModuleFP)
```

```@julia
non_flat_locus(M::ModuleFP)
```

## Depth

```@julia
koszul_homology(V::Vector, M:ModuleFP, i::Int)
```

```@julia
depth(M::ModuleFP, I::Ideal)
```

