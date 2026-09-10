# Injective Resolutions
Let $k[Q]$ be a monoid algebra for $Q\subset \mathbb{Z}^d$ and, let $M$ be a finitely 
generated $\mathbb{Z}^d$-graded module over $k[Q]$. 
An *injective resolution* is an exact sequence

$0 \to M \xhookrightarrow{\epsilon} J^0 \xrightarrow{d^0} J^1 \xrightarrow{d^1} \dots \xrightarrow{d^{i-1}} J^i \xrightarrow{d^i} \cdots.$

The maps $d^j$ are given by monomial matrices. The function [`injective_resolution`](@ref) computes 
an injective resolution up to some given cohomological degree. 
This is an implementation of the algorithms in [HM05](@cite).

!!! note
    We require that the monoid algebra $k[Q]$ is normal. 

```@docs
injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
injective_resolution(I::Oscar.InjectiveResolutions.MonoidAlgebraIdeal, i::Int)
```

### Data associated to injective modules
*Injective modules* in the $\mathbb{Z}^d$-graded setting have the form

$J = \bigoplus_{i=1}^r k\{a_i + F_i - Q\},$

where $k\{a_i + F_i - Q\}$ are so-called *indecomposable injective modules* (see, e.g., Chapter 11 of [MS05](@cite)). 

Let `J` be a $\mathbb{Z}^d$-graded injective module

$J = \bigoplus_{i=1}^r J_i = \bigoplus_{i=1}^r k\{a_i + F_i - Q\}.$

Then
- `indecomposable_injectives(J)` returns $J_1,\dots,J_r$,
- `Q_graded_part(J)` returns $J_Q = \bigoplus_{i=1}^r k\{a_i + F_i - Q\}_Q = \bigoplus_{i=1}^r k[Q]/W_i$ for some irreducible ideals $W_1,\dots,W_r$ as a finitely generated module, and
- `monoid_algebra(J)` returns the monoid algebra $k[Q]$.

### Data associated to injective resolutions
Let `I = injective_resolution(M, i)` be an injective resolution (up to cohomological degree `i`)

$0 \to M \xhookrightarrow{\epsilon} I^0 \xrightarrow{d^0} I^1 \xrightarrow{d^1} \cdots \xrightarrow{d^{i-1}} I^i$

of a finitely generated $\mathbb{Z}^d$-graded module $M$. Then
- `injective_modules(I)` returns the injective modules $I^0,I^1,\dots,I^i$,
- `cochain_maps(I)` returns the matrices of the cochain maps $d^0,d^1,\dots,d^{i-1}$,
- `monomial_matrix(j, I)` returns the differential $d^j$ as a monomial matrix, that is, together with the labels of its rows and columns,
- `degree_shift(I)` returns the degree $\alpha \in \mathbb{Z}^d$ by which $M$ was shifted internally, and
- `Q_graded_part(I)` returns the irreducible resolution of $M(-\alpha)$ from which `I` was computed.

```@docs
monomial_matrix(i::Int, res::IrrRes)
```

### Injective hulls
```@docs
injective_hull(M::SubquoModule{<:MonoidAlgebraElem})
```

### Bass numbers and minimality
```@docs
graded_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ, i::Int)
degrees_of_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
is_minimal(res::InjRes)
```

## Irreducible Resolutions
Let $M$ be a finitely generated $\mathbb{Z}^d$-graded module. An *irreducible resolution* of $M$ is an exact sequence

$0 \to M \xhookrightarrow{\epsilon} \overline{W}^0 \xrightarrow{d^0} \overline{W}^1 \xrightarrow{d^1} \cdots \xrightarrow{d^{r-1}} \overline{W}^r \rightarrow \cdots,$

where

$\overline{W}^i = \bigoplus_{j=1}^{n_i} \overline{W_{i_j}} = \bigoplus_{j=1}^{n_i} k[Q]/W_{i_j}$

for irreducible ideals $W_{i_j}$. The $k[Q]$-modules $\overline{W}^i$ are called *irreducible sums*. 

Every finitely generated $Q$-graded module has a finite minimal irreducible resolution, i.e., it is finite in length and the components are finite direct sums. It is unique up to isomorphism and obtained as the $Q$-graded part of a minimal injective resolution. For more details see, e.g., Chapter 11 of [MS05](@cite).

```@docs
irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Union{Int,Nothing})
```

### Data associated to irreducible resolutions
Let `I = irreducible_resolution(M)` be an irreducible resolution

$0 \to M \xhookrightarrow{\epsilon} \overline{W}^0 \xrightarrow{d^0} \overline{W}^1 \xrightarrow{d^1} \cdots \xrightarrow{d^{r-1}} \overline{W}^r$

of a $\mathbb{Z}^d$-graded module $M$. Then

- `irreducible_sums(I)` returns the irreducible sums $\overline{W}^0, \dots, \overline{W}^r$,
- `cochain_maps(I)` returns the cochain maps $d^0,\dots,d^{r-1}$ as module homomorphisms,
- `cochain_complex(I)` returns the complex as a `ComplexOfMorphisms{OFPModule}`,
- `is_exact(I)` checks exactness of that complex, and
- `monomial_matrix(j, I)` returns the differential $d^j$ as a monomial matrix.

```@docs
irreducible_hull(Mi::SubquoModule{<:MonoidAlgebraElem}, j)
```
