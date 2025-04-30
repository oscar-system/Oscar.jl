```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# Operations on Module Maps

If module homomorphisms `a` and `b` with `codomain(a) === domain(b)` are given,
then `compose(a, b)` refers to the composition `b` $\circ$ `a`. If an isomorphism of modules
`a` is given, then `inv(a)` refers to its inverse.

```@docs
hom_product(M::SparseFPModule, N::SparseFPModule, A::Matrix{<:SparseFPModuleHom{<:SparseFPModule, <:SparseFPModule, Nothing}})
```

```@docs
hom_tensor(M::SparseFPModule, N::SparseFPModule, V::Vector{<:SparseFPModuleHom})
```

```@docs
lift_homomorphism_contravariant(Hom_MP::SparseFPModule, Hom_NP::SparseFPModule, phi:: SparseFPModuleHom)
```

```@docs
lift_homomorphism_covariant(Hom_PM::SparseFPModule, Hom_PN::SparseFPModule, phi:: SparseFPModuleHom)
```


