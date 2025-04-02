```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
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
abelian_closure(::QQField)
```

Given the abelian closure, the generator can be recovered as follows:

```@docs
gen(::QQAbField{AbsSimpleNumField})
atlas_irrationality
atlas_description
```

## Natural embedding

Oscar assumes a natural embedding of the field `K` produced by
`K, z = abelian_closure(QQ)` into `F = algebraic_closure(QQ)`,
which is given by mapping the `n`-th root of unity returned by `z(n)`
to `root_of_unity(F, n)`.
Both roots of unity correspond to the complex number $\exp(2 \pi i / n)$.

We can convert elements of `K` to elements of `F` as follows.

```jldoctest naturalembedding
julia> K, z = abelian_closure(QQ)
(Abelian closure of rational field, Generator of abelian closure of rational field)

julia> F = algebraic_closure(QQ)
Algebraic closure of rational field

julia> x = z(5)
zeta(5)

julia> y = F(x)
Root 0.309017 + 0.951057*im of x^4 + x^3 + x^2 + x + 1

julia> y^5
Root 1.00000 of x - 1
```

Real elements of `K` can be compared with `<` and `>`.

```jldoctest naturalembedding
julia> a = x + x^4
-zeta(5)^3 - zeta(5)^2 - 1

julia> a > 0
true
```

## Printing

The n-th primitive root of the abelian closure of $\mathbf{Q}$
will by default be printed as `zeta(n)`.
The printing can be manipulated using the following functions:

```@docs
gen(::QQAbField, ::String)
set_variable!(::QQAbField, ::String)
get_variable(::QQAbField)
```

### Examples

```
julia> K, z = abelian_closure(QQ);

julia> z(4)
zeta(4)

julia> ζ = gen(K, "ζ")
Generator of abelian closure of QQ

julia> ζ(5) + ζ(3)
ζ(15)^5 + ζ(15)^3

julia> z(4)
"ζ"

julia> set_variable!(K, "zeta");

julia> z(4)
zeta(4)
```

## Reduction to characteristic ``p``

```@docs
reduce(val::QQAbFieldElem, F::FinField)
```
