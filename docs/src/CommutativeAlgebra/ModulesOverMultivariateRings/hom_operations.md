```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["hom_operations.md"]
```

# Operations on Module Maps

If module homomorphisms `a` and `b` with `codomain(a) === domain(b)` are given,
then `compose(a, b)` refers to the composition `b` $\circ$ `a`. If an isomorphism of modules
`a` is given, then `inv(a)` refers to its inverse.

```@docs
hom_product(M::ModuleFP, N::ModuleFP, A::Matrix{<: ModuleFPHom})
```

```@docs
hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{ <: ModuleFPHom})
```

```@docs
lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, phi:: ModuleFPHom)
```

```@docs
lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi:: ModuleFPHom)
```


