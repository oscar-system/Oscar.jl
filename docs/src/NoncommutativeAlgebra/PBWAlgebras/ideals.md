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

The OSCAR type for ideals in PBW-algebras is of parametrized form `PBWAlgIdeal{D, T, S}`,
where `T` is the element type of the field over which the PBW-algebra is defined (the
type `S` is added for internal use) and `D` encodes the direction
left, right, or two-sided.

## Constructors

```@docs
left_ideal(g::Vector{<:PBWAlgElem})
two_sided_ideal(g::Vector{<:PBWAlgElem})
right_ideal(g::Vector{<:PBWAlgElem})
```
