```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ideals.md"]
```

# Ideals in PBW-algebras

## Types

The OSCAR type for ideals in PBW-algebras is of parametrized form `PBWAlgIdeal{T, S}`,
where `T` is the element type of the field over which the PBW-algebra is defined (the
type `S` is added for internal use).

## Constructors

```@docs
ideal(g::Vector{<:PBWAlgElem}; two_sided=false)
```
