```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Local Cohomology
Let $M$ be a finitely generated $\mathbb{Z}^d$-graded module over a monoid algebra $k[Q]$. Further, let $J\subseteq k[Q]$ be an ideal. The $i$-th local cohomology module of $M$, denoted $H^i_J(M)$, is obtained as follows:

Let

$I^\bullet \colon 0 \to M \xrightarrow{\epsilon} I^0 \xrightarrow{d^0} I^1 \xrightarrow{d^1}\cdots \xrightarrow{d^{i-1}} I^i \xrightarrow{d^i} \cdots$

be an injective resolution of $M$. Applying the left exact functor $\Gamma_J$, that maps a $\mathbb{Z}^d$-graded module $N$ to the submodule

$\Gamma_I(N) = \{n \in \mid \exists n \in \mathbb{N} \colon n\cdot J^n = 0\},$

to $I^\bullet$ we obtain the complex

$\Gamma_J(I^\bullet) \colon 0 \to \Gamma_J(I^0) \xrightarrow{d^0} \Gamma_J(I^1) \xrightarrow{d^1}\cdots \xrightarrow{d^{i-1}} \Gamma_J(I^i) \xrightarrow{d^i} \cdots.$

The *$i$-th local cohomology module of $M$ supported on $J$* is the $i$-th cohomology module of $\Gamma_J(I^\bullet)$.

!!! note
    We require that the monoid algebra $k[Q]$ is normal. 

## Cohomological degree zero
The zeroth local cohomology module of $M$ supported by $J$ is

$H^0_J(M) = \Gamma_J(M) = \{m\in M \mid \exists n \in \mathbb{N} \colon m\cdot J^n = 0\}.$

The function `zeroth_local_cohomology` returns this submodule. 

```@docs
zeroth_local_cohomology(M::SubquoModule{<:MonoidAlgebraElem}, I::MonoidAlgebraIdeal)
```

## Sector Partition of Local Cohomology Module
The local cohomology module $H^i_J(M)$ are not finitely generated generated for $i>0$. However, sector partitions give us a finite data structure for local cohomology modules.

A *sector partition* $\mathcal{S}$ of $H^i_J(M)$ consists of

- a finite partition $\mathbb{Z}^d = \sqcup_{S\in \mathcal{S}} S$ into *sectors*
- finite dimensional $k$-vector spaces $H_S$ for each sector $S\in \mathcal{S}$, and,
- maps between these vector spaces.

Now given $\alpha \in \mathbb{Z}^d$,

$H^i_J(M)_\alpha \cong k^{\dim(H_S)} \text{ for } \alpha \in S.$

For more details on sector partitions see, e.g., Chapter 13 of [MS05](@cite).

The function [local_cohomology](@ref) computes a sector partition of $H^i_I(M)$. For performance, multiple local cohomology modules $H^1_I(M),\dots,H^i_I(M)$ should be computed using the function `local_cohomology_all`.

```@docs
    local_cohomology(M::SubquoModule{T}, I::MonoidAlgebraIdeal, i::Integer) where {T<:MonoidAlgebraElem}
    local_cohomology_all(M::SubquoModule{T}, I::MonoidAlgebraIdeal, i::Integer) where {T<:MonoidAlgebraElem}
```

### Data asssociated to Sector Partitions
Let `H = local_cohomology(M,I,i)` be a sector partition of the local cohomology module $H^i_I(M)$. Then

- `H.M` refers to $M$,
- `H.I` refers to $I$,
- `H.i` refers to `i`,
- `H.sectors` refers to the finite partition of $\mathbb{Z}^d$ into sectors as polyhedron, and,
- `H.maps` refers to the maps between the finite dimensional vector spaces. 

Each sector `S` of a sector partition consists of

- the finite dimensional $k$-vector space $H_S$ = `S.H`,
- the sector as a polyhedron `S.sector`. 

### Tests on Local Cohomology Modules
Test if a local cohomology module vanishes. 
```@docs
is_zero(S::Oscar.InjectiveResolutions.SectorPartitionLC)
```
