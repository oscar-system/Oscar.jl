```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Operations on Module Maps

If module homomorphisms `a` and `b` with `codomain(a) === domain(b)` are given,
then `compose(a, b)` refers to the composition `b` $\circ$ `a`. If an isomorphism of modules
`a` is given, then `inv(a)` refers to its inverse.

```@docs
hom_product(M::OFPModule, N::OFPModule, A::Matrix{<:OFPModuleHom{<:OFPModule, <:OFPModule, Nothing}})
```

```@docs
hom_tensor(M::OFPModule, N::OFPModule, V::Vector{<:OFPModuleHom})
```

```@docs
lift_homomorphism_contravariant(Hom_MP::OFPModule, Hom_NP::OFPModule, phi:: OFPModuleHom)
```

```@docs
lift_homomorphism_covariant(Hom_PM::OFPModule, Hom_PN::OFPModule, phi:: OFPModuleHom)
```


