```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Creating Multivariate Rings

In this section, for the convenience of the reader, we recall from the chapters on
[rings](@ref rings) and [fields](@ref fields) how to create multivariate polynomial
rings and their elements, adding illustrating examples.
At the same time, we introduce and illustrate a ring type for modelling multivariate polynomial
rings with gradings.

## Types

OSCAR provides types for dense univariate and sparse multivariate polynomials. The univariate
ring types belong to the abstract type `PolyRing{T}`, their elements have abstract type
`PolyRingElem{T}`. The multivariate ring types belong to the abstract type `MPolyRing{T}`,
their elements have abstract type `MPolyRingElem{T}`. Here, `T` is the element type
of the coefficient ring of the polynomial ring.

## Constructors

The basic constructor below allows one to build multivariate polynomial rings:

```@julia
polynomial_ring(C::Ring, xs::AbstractVector{<:VarName}; cached::Bool = true)
```

Given a ring `C` and a vector `xs` of  Symbols, Strings, or Characters, return 
a tuple `R, vars`, say, which consists of a polynomial ring `R` with coefficient ring `C`
and a vector `vars` of generators (variables) which print according to the entries of `xs`.

!!! note
    Caching is used to ensure that a given ring constructed from given parameters is unique in the system. For example, there is only one ring of multivariate polynomials over  $\mathbb{Z}$ with variables printing as x, y, z.

###### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(ZZ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over ZZ, ZZMPolyRingElem[x, y, z])

julia> typeof(R)
ZZMPolyRing

julia> typeof(x)
ZZMPolyRingElem

julia> S, (a, b, c) = polynomial_ring(ZZ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over ZZ, ZZMPolyRingElem[x, y, z])

julia> T, _ = polynomial_ring(ZZ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over ZZ, ZZMPolyRingElem[x, y, z])

julia> R === S === T
true

```

```jldoctest
julia> R1, _ = polynomial_ring(ZZ, [:x, :y, :z]);

julia> R2, _ = polynomial_ring(ZZ, ["x", "y", "z"]);

julia> R3, _ = polynomial_ring(ZZ, ['x', 'y', 'z']);

julia> R1 === R2 === R3
true

```

```jldoctest
julia> R1, x = polynomial_ring(QQ, [:x])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[x])

julia> typeof(x)
Vector{QQMPolyRingElem} (alias for Array{QQMPolyRingElem, 1})

julia> R2, (x,) = polynomial_ring(QQ, [:x])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[x])

julia> typeof(x)
QQMPolyRingElem

julia> R3, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over QQ, x)

julia> typeof(x)
QQPolyRingElem

```

```jldoctest
julia> T, x = polynomial_ring(GF(3), ["x[1]", "x[2]"]);

julia> x
2-element Vector{FqMPolyRingElem}:
 x[1]
 x[2]

```

The constructor illustrated below allows for the convenient handling of variables with multi-indices:

```jldoctest
julia> R, x, y, z = polynomial_ring(QQ, :x => (1:3, 1:4), :y => 1:2, :z => (1:1, 1:1, 1:1));

julia> x
3×4 Matrix{QQMPolyRingElem}:
 x[1, 1]  x[1, 2]  x[1, 3]  x[1, 4]
 x[2, 1]  x[2, 2]  x[2, 3]  x[2, 4]
 x[3, 1]  x[3, 2]  x[3, 3]  x[3, 4]

julia> y
2-element Vector{QQMPolyRingElem}:
 y[1]
 y[2]

julia> z
1×1×1 Array{QQMPolyRingElem, 3}:
[:, :, 1] =
 z[1, 1, 1]

```

## Coefficient Rings

Gröbner and standard bases are implemented for multivariate polynomial rings over the fields and rings below:

### The field of rational numbers $\mathbb{Q}$

```jldoctest
julia> QQ
Rational field

```
### Finite fields $\mathbb{F_p}$, $p$ a prime

```jldoctest
julia> GF(3)
Prime field of characteristic 3

julia> GF(ZZ(2)^127 - 1)
Prime field of characteristic 170141183460469231731687303715884105727

```

### Finite fields $\mathbb{F}_{p^n}$ with $p^n$ elements, $p$ a prime

```jldoctest
julia> finite_field(2, 70, "a")
(Finite field of degree 70 and characteristic 2, a)

```

### Simple algebraic extensions of $\mathbb{Q}$ or $\mathbb{F}_p$

```jldoctest
julia> T, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over QQ, t)

julia> K, a = number_field(t^2 + 1, "a")
(Number field of degree 2 over QQ, a)

julia> F = GF(3)
Prime field of characteristic 3

julia> T, t = polynomial_ring(F, :t)
(Univariate polynomial ring in t over F, t)

julia> K, a = finite_field(t^2 + 1, "a")
(Finite field of degree 2 and characteristic 3, a)

```

### Purely transcendental extensions of $\mathbb{Q}$ or $\mathbb{F}_p$

```jldoctest
julia> T, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over QQ, t)

julia> QT = fraction_field(T)
Fraction field
  of univariate polynomial ring in t over QQ

julia> parent(t)
Univariate polynomial ring in t over QQ

julia> parent(1//t)
Fraction field
  of univariate polynomial ring in t over QQ

julia> T, (s, t) = polynomial_ring(GF(3), [:s, :t]);

julia> QT = fraction_field(T)
Fraction field
  of multivariate polynomial ring in 2 variables over GF(3)

```

### The ring of integers $\mathbb{Z}$

```jldoctest
julia> ZZ
Integer ring

```

## Gradings

Given a polynomial ring $R = C[x_1, \dots, x_n]$, we may endow $R$ with various gradings.
The *standard $\mathbb Z$-grading*  on $R$ is the decomposition
$R=\bigoplus_{d\in \mathbb Z} R_d=\bigoplus_{d\geq 0} R_d$ by the usual degree of polynomials.
More general $\mathbb Z$-gradings are obtained by assigning integer weights to the variables
and considering the corresponding weighted degrees. Even more generally, we may consider
multigradings: Given a finitely generated abelian group $G$, a *multigrading* on $R$ by $G$,
or a *$G$-grading*, or simply a *grading*, corresponds to a semigroup homomorphism
$\phi: \mathbb N^n \to G$: Given $\phi$, the *degree* of a monomial $x^\alpha$
is the image $\deg(x^\alpha):=\phi(\alpha)\in G$; the induced $G$-grading on $R$
is the decomposition $R = \bigoplus_{g\in G} R_g$ satisfying $R_g\cdot R_h\subset R_{g+h}$,
where $R_g$ is the free $C$-module generated by the monomials of degree $g$. This grading is determined by
assigning the *weights* $\deg(x_i)$ to the $x_i$. In other words, it is determined by  the *weight
vector* $W = (\deg(x_1), \dots, \deg(x_n))\in G^n.$

We refer to the textbooks [MS05](@cite) and [KR05](@cite) for details on multigradings. With respect to notation,
we follow the former book.

!!! note
    Given a $G$-grading on $R$, we refer to $G$ as the *grading group* of $R$. Moreover, we then say
    that $R$ is *$G$-graded*, or simply that $R$ is *graded*.
    If $R$ is a polynomial ring over a field, we say that a $G$-grading on $R$ is *positive* if $G$ is free
    and each graded part $R_g$, $g\in G$, has finite dimension. We then also say that $R$ is
    *positively graded (by $G$)*. Note that the positivity condition can be equivalently expressed by
    asking that $G$ is free and that the degree zero part consists of the constants only (see Theorem 8.6 in [MS05](@cite)).

!!! note
    Given a  `G`-grading on `R` in OSCAR, we say that `R` is *$\mathbb Z^m$-graded* if `is_free(G) && number_of_generators(G) == rank(G) == m`
    evaluates to `true`. In this case, conversion routines allow one to switch back and forth between elements
    of `G` and integer vectors of length `m`. Specifically, if `R` is *$\mathbb Z$-graded*, that is,
    `is_free(G) && number_of_generators(G) == rank(G) == 1` evaluates to `true`,  elements of `G` may be converted
    to integers and vice versa.

### Types

Multivariate rings with gradings are modeled by objects of type
`MPolyDecRing{T, S}  :< MPolyRing{T}`, with elements of type
`MPolyRingElem_dec{T, S}  :< MPolyRingElem{T}`. Here, `S` is the element type of the
multivariate ring, and  `T` is the element type of its coefficient ring as above.

!!! note
    The types `MPolyDecRing{T, S}` and `MPolyRingElem_dec{T, S}` are
    also meant to eventually model multivariate rings with filtrations
    and their elements.

The following function allows one, in particular, to distinguish between graded and filtered rings.

```@docs
is_graded(R::MPolyRing)
```

### Constructors for Graded Rings

There are two basic ways of creating multivariate rings with gradings:
While the `grade` function allows one to create a graded ring by assigning a grading to a polynomial ring already constructed,
the `graded_polynomial_ring` function is meant to create a graded polynomial ring all at once.

```@docs
grade(R::MPolyRing, W::Vector{FinGenAbGroupElem})
```

```@docs
grade(R::MPolyRing, W::Vector{<:Vector{<:IntegerUnion}})
```

```@docs
grade(R::MPolyRing, W::Vector{<:IntegerUnion})
```
```@docs
graded_polynomial_ring(C::Ring, V::Vector{String}, W)
```

## Tests on Graded Rings


```@docs
is_standard_graded(R::MPolyDecRing)
```

```@docs
is_z_graded(R::MPolyDecRing)
```

```@docs
is_zm_graded(R::MPolyDecRing)
```

```@docs
is_positively_graded(R::MPolyDecRing)
```

## Data Associated to Multivariate Rings

Given a multivariate polynomial ring `R` with coefficient ring `C`,

- `coefficient_ring(R)` refers to `C`,
- `gens(R)` to the generators (variables) of `R`,
- `number_of_generators(R)` / `ngens(R)` to the number of these generators, and
- `gen(R, i)` as well as `R[i]` to the `i`-th such generator.

###### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> coefficient_ring(R)
Rational field

julia> gens(R)
3-element Vector{QQMPolyRingElem}:
 x
 y
 z

julia> gen(R, 2)
y

julia> R[3]
z

julia> number_of_generators(R)
3

```

In the graded case, we additionally have:

```@docs
grading_group(R::MPolyDecRing)
```

```@docs
weights(R::MPolyDecRing)
```

```@docs
monomial_basis(R::MPolyDecRing, g::FinGenAbGroupElem)
```

```@docs
homogeneous_component(R::MPolyDecRing, g::FinGenAbGroupElem)
```

```@docs
forget_grading(R::MPolyDecRing)
```

## Elements of Multivariate Rings

### Constructors

One way to create elements of a multivariate  polynomial ring is
to build up polynomials from the generators (variables) of the ring using
basic arithmetic as shown below:

###### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> f = 3*x^2+y*z
3*x^2 + y*z

julia> typeof(f)
QQMPolyRingElem

julia> S, (x, y, z) = grade(R)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> g = 3*x^2+y*z
3*x^2 + y*z

julia> typeof(g)
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}

julia> g == S(f)
true

```

Alternatively, there is the following constructor:

```@julia
(R::MPolyRing{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T <: RingElem
```

Its return value is the element of  `R`  whose nonzero coefficients are specified by the elements of `c`,
with exponent vectors given by the elements of `e`.

###### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> f = 3*x^2+y*z
3*x^2 + y*z

julia> g = R(QQ.([3, 1]), [[2, 0, 0], [0, 1, 1]])
3*x^2 + y*z

julia> f == g
true

```

An often more effective way to create polynomials is to use the `MPoly` build context as indicated below:

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> B = MPolyBuildCtx(R)
Builder for an element of R

julia> for i = 1:5 push_term!(B, QQ(i), [i, i-1]) end

julia> finish(B)
5*x^5*y^4 + 4*x^4*y^3 + 3*x^3*y^2 + 2*x^2*y + x

```


### Special Elements

Given a multivariate polynomial ring `R`, `zero(R)` and `one(R)` refer to the additive and multiplicative identity of `R`, respectively.
Relevant test calls on an element `f` of `R` are  `iszero(f)` and `isone(f)`.


### Data Associated to Elements of Multivariate Rings

Given an element `f` of a multivariate polynomial ring `R` or a graded version of such a ring, 
- `parent(f)` refers to `R`, and
- `total_degree(f)` to the total degree of `f`.

!!! note
    Given a set of variables $x = \{x_1, \ldots, x_n\}$, the *total degree* of a monomial $x^\alpha=x_1^{\alpha_1}\cdots x_n^{\alpha_n}\in\text{Mon}_n(x)$
    is the sum of the $\alpha_i$. The *total degree* of a polynomial `f`  is the maximum of the total degrees of its monomials. In particular,
    the notion of total degree ignores the weights given to the variables in the graded case.

For iterators which allow one to recover the monomials  (terms, $\dots$) of `f` we refer to the
subsection [Monomials, Terms, and More](@ref monomials_terms_more) of the section on [Gröbner/Standard Bases](@ref gb_fields).

###### Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(GF(5), [:x, :y])
(Multivariate polynomial ring in 2 variables over GF(5), FqMPolyRingElem[x, y])

julia> c = map(GF(5), [1, 2, 3])
3-element Vector{FqFieldElem}:
 1
 2
 3

julia> e = [[3, 2], [1, 0], [0, 1]]
3-element Vector{Vector{Int64}}:
 [3, 2]
 [1, 0]
 [0, 1]

julia> f = R(c, e)
x^3*y^2 + 2*x + 3*y

julia> parent(f)
Multivariate polynomial ring in 2 variables x, y
  over prime field of characteristic 5

julia> total_degree(f)
5
```

Further functionality is available in the graded case:

```@docs
homogeneous_components(f::MPolyDecRingElem{T, S}) where {T, S}
```

```@docs
homogeneous_component(f::MPolyDecRingElem, g::FinGenAbGroupElem)
```

```@docs
is_homogeneous(f::MPolyDecRingElem)
```

```@docs
degree(f::MPolyDecRingElem)
```

```@docs
forget_grading(f::MPolyDecRingElem)
```

## Homomorphisms From Multivariate Rings

If $R$ is a multivariate polynomial ring, and $S$ is any ring, then a ring homomorphism
$R \to S$ is determined by specifying its restriction to the coefficient ring of $R$,
and by assigning an image to each variable of $R$.
In OSCAR, such homomorphisms are created by using the following constructor:

```@docs
hom(R::MPolyRing, S::NCRing, coeff_map, images::Vector; check::Bool = true)
```

Given a ring homomorphism `F` from `R` to `S` as above, `domain(F)` and `codomain(F)`
refer to `R` and `S`, respectively.

!!! note
    The OSCAR homomorphism type `AffAlgHom` models ring homomorphisms `R` $\to$ `S` such that
    the type of both `R` and `S`  is a subtype of `Union{MPolyRing{T}, MPolyQuoRing{U}}`, where `T <: FieldElem` and
    `U <: MPolyRingElem{T}`. Functionality for these homomorphism is discussed in the section on [affine algebras](@ref affine_algebras).

