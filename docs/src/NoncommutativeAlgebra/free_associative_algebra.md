```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Free Associative Algebras

## Two-sided ideals

### Types

The OSCAR type for two-sided ideals in a free associative algebra is
`FreeAssociativeAlgebraIdeal{T}`, where `T` is the element type of the algebra.

### Constructors

```julia
ideal(R::FreeAssociativeAlgebra, g::Vector{T}) where T <: FreeAssociativeAlgebraElem
ideal(g::Vector{T}) where T <: FreeAssociativeAlgebraElem
```

### Ideal Membership

Non-commutative polynomial rings are not Noetherian.  Hence, in general, Groebner bases do not exist.  Hence calling the functions below may not terminate.  Picking suitable term orders is difficult in the noncommutative case.  Therefore, we fix the term order to be degree reverse lexicographic.

Setting the parameter `deg_bound` to a positive value yields the truncation of the Groebner bases to a fixed degree.  Such a truncation is always finite.

```@docs
groebner_basis(I::FreeAssociativeAlgebraIdeal, deg_bound::Int=-1; protocol::Bool=false)
```

If a finite Gröbner basis exists, it solves the ideal membership problem.

```@docs
ideal_membership(a::FreeAssociativeAlgebraElem, I::FreeAssociativeAlgebraIdeal, deg_bound::Int)
```
