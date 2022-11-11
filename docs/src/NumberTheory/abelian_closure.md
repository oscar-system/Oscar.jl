```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["abelian_closure.md"]
```

# Abelian closure of the rationals

The abelian closure $\mathbf{Q}^\text{ab}$ is the maximal abelian extension of $\mathbf{Q}$
inside a fixed algebraic closure and can explicitly described as
```math
\mathbf{Q}^{\mathrm{ab}} = \mathbf{Q}(\zeta_n \mid n \in \mathbf{N}),
```
the union of all cyclotomic extensions. Here for $n \in \mathbf{N}$ we denote by $\zeta_n$ a primitive $n$-th root of unity.

## Creation of the abelian closure and elements

```@docs
abelian_closure(::FlintRationalField)
```

Given the abelian closure, the generator can be recovered as follows:

```@docs
gen(::QQAbField)
```

## Printing

The n-th primitive root of the abelian closure of will by default be printed as
`z(n)`. The printing can be manipulated using the following functions:

```@docs
gen(::QQAbField, ::String)
set_variable!(::QQAbField, ::String)
get_variable(::QQAbField)
```

### Examples

```@jldoctest
julia> K, z = abelian_closure(QQ);

julia> z(4)
z(4)

julia> ζ = gen(K, "ζ")
Generator of abelian closure of Q

julia> ζ(5) + ζ(3)
ζ(15)^5 + ζ(15)^3
```
