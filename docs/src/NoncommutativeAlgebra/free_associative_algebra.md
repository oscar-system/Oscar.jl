```@meta
CurrentModule = Oscar
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

Non-commutative polynomial rings are not Noetherian.  Hence, in general, Groebner bases do not exist.  Hence calling the functions below may not terminate.  Picking suitable term orders is difficult in the noncommutative case.  Therefore, we fix the term order to be degree reverse lexicographic.

Setting the parameter `deg_bound` to a positive value yields the truncation of the Groebner bases to a fixed degree.  Such a truncation is always finite.

```@docs
groebner_basis(I::FreeAssAlgIdeal, deg_bound::Int=-1; protocol::Bool=false)
```

If a finite Gröbner basis exists, it solves the ideal membership problem.

```@docs
ideal_membership(a::FreeAssAlgElem, I::FreeAssAlgIdeal, deg_bound::Int)
```
