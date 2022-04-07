```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["free_associative_algebra.md"]
```

# Free Associative Algebras

## Two-sided ideals

### Types

The OSCAR type for two-sided ideals in a free associative algebra is
`FreeAssAlgIdeal{T}`, where `T` is the element type of the algebra.

### Constructors

```julia
ideal(R::FreeAssAlgebra, g::Vector{T}) where T <: FreeAssAlgElem
ideal(g::Vector{T}) where T <: FreeAssAlgElem
```

### Ideal Membership

```@docs
ideal_membership(a::FreeAssAlgElem, I::FreeAssAlgIdeal, deg_bound::Int)
```

