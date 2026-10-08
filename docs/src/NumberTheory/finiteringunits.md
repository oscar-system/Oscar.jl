```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Unit group and first K-group

We have the following experimental functions to determine unit groups and the first K-group of finite rings and finite-dimensional algebras.

```@docs
unit_group(::Union{FiniteRing, Oscar.Hecke.AbstractAssociativeAlgebra})
Oscar.k1(::Union{FiniteRing, Oscar.Hecke.AbstractAssociativeAlgebra})
```
