```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Injective Resolutions
Let $k[Q]$ be a monoid algebra for $Q\subset \mathbb{Z}^d$ and, let $M$ be a finitely 
generated $\mathbb{Z}^d$-graded module over $k[Q]$. 
An *injective resolution* is an exact sequence

$0 \to M \xhookrightarrow{\epsilon} J^0 \xrightarrow{d^0} J^1 \xrightarrow{d^1} \dots \xrightarrow{d^{i-1}} J^i \xrightarrow{d^i} \cdots.$

The maps $d^j$ are given by monomial matrices. The function [injective_resolution](@ref) computes 
an injective resolution up to some given cohomological degree. 
This is an implementation of the algorithms in [HM05](@cite).

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
- `J.Q_graded_part` refers to $J_Q = \bigoplus_{i=1}^r k\{a_i + F_i - Q\}_Q = \bigoplus_{i=1}^r k[Q]/W_i$ for some irreducible ideals $W_1,\dots,W_r$,
- `J.indec_injectives` refers to $J_1,\dots,J_r$, and, 
- `J.monoid_algebra` refers to the monoid algebra $k[Q]$.

### Data associated to injective resolutions
Let `I = injective_resolution(M,i)` be an injective resolution (up to cohomological degree i)

$0 \to M \xhookrightarrow{\epsilon} I^0 \xrightarrow{d^0} I^1 \xrightarrow{d^1} \cdots \xrightarrow{d^{i-1}} I^i$

of a finitely generated $Q$-graded module $M$. Then
- `I.cochain_maps` refers to the cochain maps $d^0,d^1,d^2,\dots,d^{i-1}$,
- `I.inj_mods` refers to the injective modules $I^0,I^1,I^2,\dots,I^i$,
- `I.mod` refers to $M$,
- `I.upto` refers to the length of `I`,
- `I.shift` refers to a degree $\alpha \in \mathbb{Z}^d$, and, 
- `I.Q_graded_part` refers to the irreducible resolution of $M$ shifted by $\alpha$.

## Irreducible Resolutions
Let $M$ be a finitely generated $\mathbb{Z}^d$-graded module. An *irreducible resolution* of $M$ is an exact sequence

$0 \to M \xhookrightarrow{\epsilon} \overline{W}^0 \xrightarrow{d^0} \overline{W}^1 \xrightarrow{d^1} \cdots \xrightarrow{d^{r-1}} \overline{W}^r \rightarrow \cdots,$

where

$\overline{W}^i = \bigoplus_{j=1}^{n_i} \overline{W_{i_j}} = \bigoplus_{j=1}^{n_i} k[Q]/W_{i_j}$

for irreducible ideals $W_{i_j}$. The $k[Q]$-modules $\overline{W}^i$ are called *irreducible sums*. 

Every finitely generated $Q$-graded module has a finite minimal irreducible resolution, i.e., it is finite in length and the components are finite direct sums. It is unique up to isomorphism and obtained as the $Q$-graded part of a minimal injective resolution. For more details see, e.g., Chapter 11 of [MS05](@cite).

```@docs
irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int=0)
```

### Data associated to irreducible resolutions
Let `I = irreducible_resolution(M)` be an irreducible resolution

$0 \to M \xhookrightarrow{\epsilon} \overline{W}^0 \xrightarrow{d^0} \overline{W}^1 \xrightarrow{d^1} \cdots \xrightarrow{d^{r-1}} \overline{W}^r$

of a $\mathbb{Z}^d$-graded module $M$. Then

- `I.cochain_maps` refers to the cochain maps $d^0,\dots,d^{r-1}$,
- `I.irr_sums` refers to the irreducible sums $\overline{W}^0, \dots, \overline{W}^r$,
- `I.mod` refers to $M$, and,
- `I.cochain_complex` refer to the exact sequence as a `ComplexOfMorphisms{ModuleFP}`.